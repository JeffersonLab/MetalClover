#include "compact_solve.hpp"
#include "production_bridge.hpp"
#include "three_level.hpp"
#include <filesystem>
#include <iomanip>
#include <sys/resource.h>

namespace {
  template <class T> void release(std::vector<T> &values) {
    std::vector<T>().swap(values);
  }

  class Import {
    std::ifstream stream;
    std::streampos reference;

  public:
    production::Header header;
    explicit Import(const std::string &prefix)
        : stream(prefix + ".bridge", std::ios::binary),
          header(production::read_header(stream)) {}
    std::vector<M3> links() {
      std::vector<M3> links(size_t(header.volume()) * 4);
      std::array<chroma_bridge::Complex, 9 * 256> buffer;
      for (size_t first = 0; first < links.size(); first += 256) {
        size_t n = std::min(size_t(256), links.size() - first);
        chroma_bridge::read(stream, buffer.data(), n * 9);
        for (size_t i = 0; i < n; ++i)
          for (int k = 0; k < 9; ++k) {
            auto z = buffer[i * 9 + k];
            if (!std::isfinite(z.re) || !std::isfinite(z.im))
              throw std::runtime_error("nonfinite stout link");
            links[first + i][k] = {z.re, z.im};
          }
      }
      reference = stream.tellg();
      return links;
    }
    double compare(const std::vector<D> &values, double sign = 1) {
      stream.clear();
      stream.seekg(reference);
      std::array<chroma_bridge::Complex, 4096> buffer;
      double error = 0, norm = 0;
      for (size_t first = 0; first < values.size(); first += buffer.size()) {
        size_t n = std::min(buffer.size(), values.size() - first);
        chroma_bridge::read(stream, buffer.data(), n);
        for (size_t j = 0; j < n; ++j) {
          auto a = values[first + j];
          auto b = buffer[j];
          error += std::norm(CD(sign * a.re - b.re, sign * a.im - b.im));
          norm += std::norm(CD(b.re, b.im));
        }
      }
      if (stream.peek() != EOF)
        throw std::runtime_error("unexpected production bridge tail");
      double rel = std::sqrt(error / norm);
      if (!std::isfinite(rel) || rel > 1e-12)
        throw std::runtime_error(
            "production full operator differs from Chroma");
      return rel;
    }
    void validate(CompactClover &c, int workers) {
      c.rhs(1, false);
      std::vector<D> x(c.w.a.p.outputs()), y(x.size());
      for (size_t q = 0; q < x.size(); ++q) {
        auto z = chroma_bridge::source(q, chroma_bridge::probe_rhs);
        x[q] = {z.re, z.im};
      }
      apply_cpu(c.w.a, x.data(), y.data(), workers);
      double original = compare(y);
      c.w.a.b.resize(x.size());
      std::vector<double> sums(workers);
      compact_fp64::residual(c, x.data(), c.w.a.b.data(), y.data(), sums.data(),
                             workers);
      double fused = compare(y, -1);
      release(c.w.a.b);
      std::cout << "CHROMA MATCH original=" << original
                << " optimized_FP64=" << fused << std::endl;
    }
  };

  struct EigenZero {
    uint32_t first, count, spatial;
    std::vector<Cx<float>> values;
    EigenZero(const std::string &prefix, const production::Header &h) {
      std::ifstream in(prefix + ".eigen0", std::ios::binary);
      std::array<char, 8> tag;
      chroma_bridge::read(in, tag.data(), 8);
      if (tag != production::eigen_magic)
        throw std::runtime_error("wrong eigenvector interchange");
      for (auto *x : {&first, &count, &spatial})
        chroma_bridge::read(in, x, 1);
      if (first != h.first_time || count != h.time_count ||
          spatial != h.shape[0] * h.shape[1] * h.shape[2])
        throw std::runtime_error("eigenvector geometry mismatch");
      values.resize(size_t(count) * spatial * 3);
      chroma_bridge::read(in, values.data(), values.size());
      if (in.peek() != EOF)
        throw std::runtime_error("unexpected eigenvector tail");
    }
    const Cx<float> *at(uint32_t time) const {
      if (time < first || time >= first + count)
        throw std::runtime_error("missing eigenvector timeslice");
      return values.data() + size_t(time - first) * spatial * 3;
    }
    void source(CompactClover &c, CompactOuter &outer, uint32_t time,
                uint32_t spin, int workers) const {
      auto &b = c.w.a.b;
      std::fill(b.begin(), b.end(), D{});
      auto *eigen = at(time);
      for (uint32_t site = 0; site < spatial; ++site)
        for (uint32_t color = 0; color < 3; ++color) {
          auto z = eigen[size_t(site) * 3 + color];
          b[(size_t(time) * spatial + site) * 12 + color * 4 + spin] = {z.re,
                                                                        z.im};
        }
      compact_rhs_norms(c.w.a, outer, workers);
    }
    void project(const std::vector<D> &x, uint32_t source_spin,
                 std::ostream &out) const {
      for (uint32_t time = first; time < first + count; ++time) {
        auto *eigen = at(time);
        for (uint32_t spin = 0; spin < 4; ++spin) {
          CD sum = 0;
          for (uint32_t site = 0; site < spatial; ++site)
            for (uint32_t color = 0; color < 3; ++color) {
              auto v = eigen[size_t(site) * 3 + color];
              auto z =
                  x[(size_t(time) * spatial + site) * 12 + color * 4 + spin];
              sum += std::conj(CD(v.re, v.im)) * CD(z.re, z.im);
            }
          out << source_spin << ',' << spin << ',' << time << ',' << sum.real()
              << ',' << sum.imag() << '\n';
        }
      }
      out.flush();
    }
  };

  void check_metal(CompactClover &c, MetalRun &metal, int workers) {
    c.rhs(1, false);
    std::vector<Cx<float>> x(c.parity.oe.p.outputs()), y(x.size());
    for (size_t q = 0; q < x.size(); ++q) {
      auto z = chroma_bridge::source(q, 4321);
      x[q] = {float(z.re), float(z.im)};
    }
    CPUCloverSchur<float, CompactClover> cpu(c, workers);
    cpu.apply(x.data(), y.data());
    GPUCompactSchur gpu(metal, c);
    auto in = [metal.device newBufferWithBytes:x.data()
                                        length:x.size() * 8
                                       options:MTLResourceStorageModeShared];
    auto out = [metal.device newBufferWithLength:x.size() * 8
                                         options:MTLResourceStorageModeShared];
    auto cmd = [metal.queue commandBuffer];
    auto enc = [cmd computeCommandEncoder];
    gpu.encode(enc, in, out);
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
    if (cmd.status == MTLCommandBufferStatusError)
      throw std::runtime_error("production device comparison failed");
    auto *got = static_cast<Cx<float> *>([out contents]);
    double error = 0, norm = 0;
    for (size_t q = 0; q < y.size(); ++q) {
      error += std::norm(CD(got[q].re - y[q].re, got[q].im - y[q].im));
      norm += std::norm(CD(y[q].re, y[q].im));
    }
    double relative = std::sqrt(error / norm);
    if (!std::isfinite(relative) || relative > 3e-6)
      throw std::runtime_error("production CPU/Metal Schur mismatch");
    std::cout << "METAL SCHUR MATCH relative=" << relative << std::endl;
  }

  void check_saved(const CompactClover &c, CompactOuter &outer, int workers) {
    compact_fp64::residual(c, outer.x.data(), c.w.a.b.data(),
                           outer.residual.data(), outer.sums.data(), workers);
    double n = 0;
    for (double v : outer.sums)
      n += v;
    double relative = std::sqrt(n) / c.w.a.bnorm[0];
    if (!std::isfinite(relative) || relative > 1e-10)
      throw std::runtime_error("saved solution fails full FP64 residual");
    std::cout << "SAVED full_relative=" << relative << std::endl;
  }
} // namespace

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      if (argc < 3)
        throw std::runtime_error("Usage: production_solve BRIDGE_PREFIX "
                                 "OUTPUT_PREFIX [--validate-only] [--resume] "
                                 "[--sweeps N] [--initial-prefix PREFIX]");
      std::string prefix = argv[1], output = argv[2], initial_prefix;
      bool validate = false, resume = false;
      int sweeps = 8;
      for (int i = 3; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--validate-only")
          validate = true;
        else if (flag == "--resume")
          resume = true;
        else if (flag == "--sweeps" && i + 1 < argc)
          sweeps = std::stoi(argv[++i]);
        else if (flag == "--initial-prefix" && i + 1 < argc)
          initial_prefix = argv[++i];
        else
          throw std::runtime_error("invalid production solver option");
      }
      if (sweeps < 1 || sweeps > 64)
        throw std::runtime_error("invalid filter sweep count");
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
      const int workers = 8;
      Import input(prefix);
      auto h = input.header;
      std::string source;
      for (auto file :
           {"kernels.metal", "bicgstab.metal", "clover.metal",
            "multigrid.metal", "sparse_coarse.metal", "compact_clover.metal"}) {
        std::ifstream f(file);
        source.append((std::istreambuf_iterator<char>(f)), {});
        source += '\n';
      }
      MetalRun metal(source);
      std::cout << std::setprecision(13) << "GPU "
                << [[metal.device name] UTF8String] << " shape=" << h.shape[0]
                << 'x' << h.shape[1] << 'x' << h.shape[2] << 'x' << h.shape[3]
                << std::endl;
      CompactClover c(h.shape, h.mass, workers, .25, input.links(), h.xi0, h.nu,
                      h.clover_r, h.clover_t);
      release(c.evenA32);
      std::cout << "CLOVER setup_s=" << c.setup_ms / 1000
                << " min_pivot=" << c.min_pivot << std::endl;
      input.validate(c, workers);
      check_metal(c, metal, workers);
      if (validate) {
        std::cout << "VALIDATION PASS" << std::endl;
        return 0;
      }

      LatticeShape blocks{4, 4, 4, 8}, coarse_shape;
      for (int mu = 0; mu < 4; ++mu) {
        if (h.shape[mu] % blocks[mu])
          throw std::runtime_error("production block must divide geometry");
        coarse_shape[mu] = h.shape[mu] / blocks[mu];
      }
      auto setup = Clock::now();
      const uint32_t vectors = 12, filter_rhs = 2;
      c.rhs(filter_rhs, false);
      Transfer first(c.parity.oe.p.rows, 12, vectors, h.shape, blocks, true,
                     workers, &c.parity.odd);
      {
        CPUCloverSchur<double, CompactClover> schur(c, workers);
        first.build([&](const D *x, D *y) { schur.apply(x, y); },
                    [&](const D *x, D *y) {
                      local_blocks(c.inverseOdd32, x, y, first.rows, filter_rhs,
                                   workers);
                    },
                    sweeps, .65, filter_rhs);
      }
      std::cout << "P1 built coarse=" << first.coarse << std::endl;
      BlockMatrix a1(2 * vectors, first.aggregates, workers);
      assemble_schur(a1, c, first);
      release(c.inverseEven);
      std::cout << "A1 blocks=" << a1.col.size() << std::endl;
      Transfer second(a1.rows, a1.bs, vectors, coarse_shape, hypercube(2),
                      false, workers);
      second.build([&](const D *x, D *y) { a1.apply(x, y, vectors); },
                   [&](const D *x, D *y) { a1.apply(x, y, vectors, true); },
                   sweeps);
      BlockMatrix a2(2 * vectors, second.aggregates, workers);
      project_matrix(a2, a1, second);
      first.release64();
      second.release64();
      a1.release64();
      a2.release64();
      std::cout << "HIERARCHY built seconds=" << seconds(setup)
                << " coarse1=" << a1.n << " coarse2=" << a2.n << std::endl;
      c.rhs(1, false);
      c.w.a.b.resize(c.w.a.p.outputs());
      c.w.a.bnorm.resize(1);
      CompactOuter outer(c.w.a, workers);
      outer.known_solution = false;
      EigenZero eigen(prefix, h);
      GPUCompactSchur op(metal, c);
      release(c.links32);
      release(c.oddA32);
      release(c.inverseEven32);
      GPUCycle pre(first, metal, 1, c.inverseOdd32);
      release(first.basis32);
      release(c.inverseOdd32);
      pre.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                        id<MTLBuffer> y) { op.encode(e, x, y); };
      GPUBlockCoarseSolve<BlockMatrix> matrix1(a1, metal, 1, 1, .1),
          bottom(a2, metal, 1, 16, .1);
      GPUBlockOperator op1(matrix1);
      GPUCycle pre1(second, metal, 1, a1.di32);
      release(second.basis32);
      pre1.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                         id<MTLBuffer> y) { op1.encode(e, x, y); };
      pre1.coarse = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                        id<MTLBuffer> y) { bottom.encode(e, x, y); };
      pre1.check = [&]() { bottom.check(); };
      GPUGenericFGMRES<GPUBlockOperator, GPUCycle> fg1(a1.n, 1, op1, pre1, 8);
      fg1.allow_inexact = true;
      pre.blocking_coarse = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
        fg1.run_device(x, .1, 16);
        pre.copy(fg1.delta, y, a1.n);
      };
      release(a1.a32);
      release(a1.di32);
      release(a2.a32);
      release(a2.di32);
      GPUGenericFGMRES<GPUCompactSchur, GPUCycle> fine(first.rows * 12, 1, op,
                                                       pre, 8);
      fine.blocking_pre = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
        pre.apply_device(x, y);
      };
      fine.progress = [](int k, double r) {
        std::cout << "INNER k=" << k << " relative_FP32=" << r << std::endl;
      };
      CompactDeviceCorrection<decltype(fine)> correction(c, op, fine, workers);
      std::cout << "METAL allocated_GiB="
                << double(metal.device.currentAllocatedSize) / 1073741824.
                << std::endl;
      if (!resume && (std::filesystem::exists(output + ".csv") ||
                      std::filesystem::exists(output + ".trace.csv")))
        throw std::runtime_error("production output already exists");
      std::ofstream projections(output + ".csv");
      projections << std::setprecision(17)
                  << "spin_source,spin_sink,time,real,imag\n";
      std::ofstream traces(output + ".trace.csv",
                           resume ? std::ios::app : std::ios::out);
      if (!resume)
        traces << "spin_source,outer,krylov,full_relative\n";
      for (uint32_t spin = 0; spin < 4; ++spin) {
        eigen.source(c, outer, h.source_time, spin, workers);
        std::string path = output + ".spin" + std::to_string(spin) + ".bin";
        if (std::filesystem::exists(path)) {
          if (!resume)
            throw std::runtime_error("solution output exists");
          std::ifstream old(path, std::ios::binary);
          chroma_bridge::read(old, outer.x.data(), outer.x.size());
          if (old.peek() != EOF)
            throw std::runtime_error("saved solution has wrong size");
          check_saved(c, outer, workers);
        } else {
          outer.observe = [&, spin](int step, double norm, int krylov) {
            traces << std::setprecision(17) << spin << ',' << step << ','
                   << krylov << ',' << norm << '\n';
            traces.flush();
            std::cout << "OUTER spin=" << spin << " step=" << step
                      << " K=" << krylov << " full_relative=" << norm
                      << std::endl;
          };
          fine.fine_calls = pre.products = 0;
          bottom.reset();
          bool initial = false;
          if (!initial_prefix.empty()) {
            std::string guess =
                initial_prefix + ".spin" + std::to_string(spin) + ".bin";
            if (std::filesystem::exists(guess)) {
              std::ifstream old(guess, std::ios::binary);
              chroma_bridge::read(old, outer.x.data(), outer.x.size());
              if (old.peek() != EOF)
                throw std::runtime_error("initial solution has wrong size");
              initial = true;
              std::cout << "INITIAL guess spin=" << spin << " path=" << guess
                        << std::endl;
            }
          }
          auto result =
              compact_solve(c, correction, outer, workers, 8192, initial);
          std::ofstream saved(path + ".tmp", std::ios::binary);
          chroma_bridge::write(saved, outer.x.data(), outer.x.size());
          saved.close();
          if (!saved)
            throw std::runtime_error("failed to save converged solution");
          std::filesystem::rename(path + ".tmp", path);
          std::cout << "SOLVED spin=" << spin
                    << " seconds=" << result.total_ms / 1000
                    << " residual_s=" << result.residual_ms / 1000
                    << " K=" << result.inner_max_sum
                    << " full_relative=" << result.max_residual << std::endl;
        }
        eigen.project(outer.x, spin, projections);
      }
      if (!projections || !traces)
        throw std::runtime_error("failed to write production results");
      std::cout << "PASS four v0 spin sources, full FP64 residuals <= 1e-10"
                << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
