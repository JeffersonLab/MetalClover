#include "mg_real_tune.hpp"
#include "compact_solve.hpp"
#include "production_bridge.hpp"

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>

using namespace MetalCloverEnv;

namespace {
  constexpr double target = 1e-12;
  constexpr int workers = 8;
  struct CountedSchur : GPUCompactSchur {
    using GPUCompactSchur::GPUCompactSchur;
    uint64_t applications = 0;
    void encode(id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                id<MTLBuffer> y) {
      ++applications;
      GPUCompactSchur::encode(e, x, y);
    }
  };
  template <class T> void release(std::vector<T> &v) {
    std::vector<T>().swap(v);
  }
  void source(CompactClover &c, CompactOuter &w, int rhs) {
    double scale = 1 / std::sqrt(double(c.w.a.p.outputs()));
    chunks(c.w.a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t q = first; q < last; ++q) {
        auto z = chroma_bridge::source(q, 20261107 + rhs);
        c.w.a.b[q] = {scale * z.re, scale * z.im};
      }
    });
    compact_rhs_norms(c.w.a, w, workers);
  }
  struct Input {
    std::ifstream file;
    production::Header h;
    std::streampos probe;
    explicit Input(const std::string &prefix)
        : file(prefix + ".bridge", std::ios::binary),
          h(production::read_header(file)) {}
    std::vector<M3> links() {
      std::vector<M3> out(size_t(h.volume()) * 4);
      std::array<chroma_bridge::Complex, 9> row;
      for (auto &m : out) {
        chroma_bridge::read(file, row.data(), 9);
        for (int k = 0; k < 9; ++k)
          m[k] = {row[k].re, row[k].im};
      }
      probe = file.tellg();
      return out;
    }
    double compare(const std::vector<D> &x) {
      file.clear();
      file.seekg(probe);
      std::vector<D> reference(x.size());
      chroma_bridge::read(file, reference.data(), reference.size());
      double e = 0, n = 0;
      for (size_t j = 0; j < x.size(); ++j) {
        e +=
            std::norm(CD(x[j].re - reference[j].re, x[j].im - reference[j].im));
        n += reference[j].re * reference[j].re +
             reference[j].im * reference[j].im;
      }
      double rel = std::sqrt(e / n);
      if (!std::isfinite(rel) || rel > 1e-13)
        throw std::runtime_error("standalone full operator mismatch");
      return rel;
    }
  };
  double full_residual(const CompactClover &c, CompactOuter &w) {
    compact_fp64::residual(c, w.x.data(), c.w.a.b.data(), w.residual.data(),
                           w.sums.data(), workers);
    double n = 0;
    for (double v : w.sums)
      n += v;
    return std::sqrt(n) / c.w.a.bnorm[0];
  }
  void check_schur(CompactClover &c, GPUCompactSchur &op) {
    std::vector<Cx<float>> x(c.parity.oe.p.outputs()), reference(x.size());
    for (size_t j = 0; j < x.size(); ++j) {
      auto z = chroma_bridge::source(j, 20261108);
      x[j] = {float(z.re), float(z.im)};
    }
    CPUCloverSchur<float, CompactClover> cpu(c, workers);
    cpu.apply(x.data(), reference.data());
    auto input = [op.m.device newBufferWithBytes:x.data()
                                          length:x.size() * 8
                                         options:MTLResourceStorageModeShared];
    auto output =
        [op.m.device newBufferWithLength:x.size() * 8
                                 options:MTLResourceStorageModeShared];
    auto cmd = [op.m.queue commandBuffer];
    auto enc = [cmd computeCommandEncoder];
    op.encode(enc, input, output);
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
    if (cmd.status == MTLCommandBufferStatusError)
      throw std::runtime_error("Schur probe GPU error");
    auto *values = static_cast<Cx<float> *>([output contents]);
    double e = 0, n = 0;
    for (size_t j = 0; j < x.size(); ++j) {
      e += std::norm(
          CD(values[j].re - reference[j].re, values[j].im - reference[j].im));
      n +=
          reference[j].re * reference[j].re + reference[j].im * reference[j].im;
    }
    double rel = std::sqrt(e / n);
    if (!std::isfinite(rel) || rel > 3e-6)
      throw std::runtime_error("CPU/Metal Schur convention mismatch");
    std::cout << "METAL CPU_FP32 Schur_match=" << rel << std::endl;
  }

  struct Options {
    std::string bridge, output, cycles = "1:1:1";
    LatticeShape blocks{3, 3, 3, 8};
    int vectors = 8, second_vectors = 8, sweeps = 128, restart = 12,
        coarse_limit = 32, bottom_limit = 4, rhs = 2, screen_outers = 0,
        screen_limit = 2048, finish = 1, bicg = 0;
    double filter_omega = .65, coarse_tol = .1, bottom_tol = .1;
  };
  LatticeShape shape_arg(const std::string &text) {
    LatticeShape a{};
    std::stringstream s(text);
    for (int i = 0; i < 4; ++i) {
      char c;
      if (i && (!(s >> c) || c != ','))
        throw std::runtime_error("four comma-separated extents required");
      if (!(s >> a[i]) || !a[i])
        throw std::runtime_error("invalid extent");
    }
    char c;
    if (s >> c)
      throw std::runtime_error("extent suffix");
    return a;
  }
  Options options(int argc, char **argv) {
    if (argc < 3 || (argc - 3) % 2)
      throw std::runtime_error(
          "Usage: mg_real_tune BRIDGE_PREFIX OUTPUT_PREFIX [--name value] ...");
    Options o;
    o.bridge = argv[1];
    o.output = argv[2];
    for (int i = 3; i < argc; i += 2) {
      std::string k = argv[i], v = argv[i + 1];
      if (k == "--blocks")
        o.blocks = shape_arg(v);
      else if (k == "--cycles")
        o.cycles = v;
      else if (k == "--vectors")
        o.vectors = std::stoi(v);
      else if (k == "--second-vectors")
        o.second_vectors = std::stoi(v);
      else if (k == "--filter-sweeps")
        o.sweeps = std::stoi(v);
      else if (k == "--filter-omega")
        o.filter_omega = std::stod(v);
      else if (k == "--restart")
        o.restart = std::stoi(v);
      else if (k == "--coarse-limit")
        o.coarse_limit = std::stoi(v);
      else if (k == "--bottom-limit")
        o.bottom_limit = std::stoi(v);
      else if (k == "--coarse-tol")
        o.coarse_tol = std::stod(v);
      else if (k == "--bottom-tol")
        o.bottom_tol = std::stod(v);
      else if (k == "--rhs")
        o.rhs = std::stoi(v);
      else if (k == "--screen-outers")
        o.screen_outers = std::stoi(v);
      else if (k == "--screen-limit")
        o.screen_limit = std::stoi(v);
      else if (k == "--finish")
        o.finish = std::stoi(v);
      else if (k == "--bicg")
        o.bicg = std::stoi(v);
      else
        throw std::runtime_error("unknown option " + k);
    }
    if (o.vectors < 1 || o.vectors > 32 || o.second_vectors < 1 ||
        o.second_vectors > 32 || o.sweeps < 0 || o.sweeps > 2048 ||
        o.restart < 2 || o.restart > 24 || o.rhs < 1 || o.rhs > 192 ||
        o.screen_outers < 0 || o.screen_outers > 3 || o.screen_limit < 1 ||
        o.coarse_limit < 1 || o.bottom_limit < 1 ||
        !(o.filter_omega > 0 && o.filter_omega <= 2) ||
        !(o.coarse_tol > 0 && o.coarse_tol < 1) ||
        !(o.bottom_tol > 0 && o.bottom_tol < 1))
      throw std::runtime_error("option out of range");
    return o;
  }
  struct Result {
    double elapsed = 0, rho = 1, previous = 1, last_seconds = 0, residual_s = 0;
    uint64_t applications = 0;
    int krylov = 0, outers = 0;
  };
  template <class Correction>
  Result refine(const CompactClover &c, CompactOuter &w, Correction &correction,
                const std::string &method, const std::string &output,
                std::ostream &trace, std::ostream &summary, int rhs,
                CountedSchur &op, int max_outer, int limit, bool save) {
    std::fill(w.x.begin(), w.x.end(), D{});
    w.residual = c.w.a.b;
    const uint64_t before = op.applications;
    auto start = Clock::now();
    Result r;
    for (int outer = 1; outer <= max_outer; ++outer) {
      auto step = Clock::now();
      const uint64_t step_before = op.applications;
      auto info = correction.run(w.residual, .002, limit);
      r.krylov += info.max_iterations;
      auto *dx = correction.correction();
      chunks(c.w.a.p.rows * 12, workers,
             [&](int, uint32_t first, uint32_t last) {
               for (size_t j = first; j < last; ++j) {
                 w.x[j].re += double(dx[j].re);
                 w.x[j].im += double(dx[j].im);
               }
             });
      auto check = Clock::now();
      r.previous = r.rho;
      r.rho = full_residual(c, w);
      r.residual_s += seconds(check);
      r.outers = outer;
      r.last_seconds = seconds(step);
      if (!std::isfinite(r.rho))
        throw std::runtime_error("nonfinite full FP64 residual");
      trace << method << ',' << rhs << ',' << outer << ','
            << info.max_iterations << ',' << r.rho << ','
            << op.applications - step_before << ',' << r.last_seconds << '\n';
      trace.flush();
      std::cout << "OUTER method=" << method << " rhs=" << rhs
                << " step=" << outer << " K=" << info.max_iterations
                << " fine_S=" << op.applications - step_before
                << " seconds=" << r.last_seconds << " full_relative=" << r.rho
                << std::endl;
      if (r.rho <= target)
        break;
    }
    r.elapsed = seconds(start);
    r.applications = op.applications - before;
    double io_seconds = 0;
    if (save) {
      if (r.rho > target)
        throw std::runtime_error("full FP64 1e-12 refinement failed");
      std::string path =
          output + "." + method + ".rhs" + std::to_string(rhs) + ".bin";
      if (std::filesystem::exists(path))
        throw std::runtime_error("solution output already exists");
      auto io = Clock::now();
      std::ofstream saved(path, std::ios::binary);
      chroma_bridge::write(saved, w.x.data(), w.x.size());
      saved.close();
      if (!saved)
        throw std::runtime_error("solution write failed");
      io_seconds = seconds(io);
    }
    summary << method << ',' << rhs << ',' << r.outers << ',' << r.krylov << ','
            << r.elapsed << ',' << r.residual_s << ',' << r.rho << ','
            << io_seconds << ',' << r.applications << '\n';
    summary.flush();
    std::cout << (save ? "SOLVED" : "SCREENED") << " method=" << method
              << " rhs=" << rhs << " outer=" << r.outers << " K=" << r.krylov
              << " fine_S=" << r.applications << " seconds=" << r.elapsed
              << " full_relative=" << r.rho << std::endl;
    return r;
  }
  struct Cycle {
    int pre, post, bottom;
    float omega;
    std::string name;
  };
  std::vector<Cycle> cycles_arg(const std::string &text, int default_bottom) {
    std::vector<Cycle> result;
    std::stringstream list(text);
    std::string item;
    while (std::getline(list, item, ',')) {
      Cycle c{};
      char sep;
      std::stringstream s(item);
      if (!(s >> c.pre >> sep) || sep != ':' || !(s >> c.post >> sep) ||
          sep != ':' || !(s >> c.omega) || c.pre < 0 || c.post < 0 ||
          c.pre > 4 || c.post > 4 || c.omega <= 0 || c.omega > 2)
        throw std::runtime_error("invalid cycle pre:post:weight");
      c.bottom = default_bottom;
      if (s >> sep) {
        if (sep != ':' || !(s >> c.bottom) || c.bottom < 1)
          throw std::runtime_error("invalid bottom cap");
      }
      c.name = "screen-" + std::to_string(result.size());
      result.push_back(c);
    }
    if (result.empty())
      throw std::runtime_error("empty cycles");
    return result;
  }
} // namespace
int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      Options o = options(argc, argv);
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
      if (std::filesystem::exists(o.output + ".metal-solves.csv"))
        throw std::runtime_error("choose a fresh output prefix");
      Input in(o.bridge);
      auto h = in.h;
      LatticeShape coarse_shape{};
      for (int mu = 0; mu < 4; ++mu) {
        if (h.shape[mu] % o.blocks[mu] || (h.shape[mu] / o.blocks[mu]) % 2)
          throw std::runtime_error(
              "blocks must divide and leave even coarse extents");
        coarse_shape[mu] = h.shape[mu] / o.blocks[mu];
      }
      CompactClover c(h.shape, h.mass, workers, .25, in.links(), h.xi0, h.nu,
                      h.clover_r, h.clover_t);
      release(c.evenA32);
      c.rhs(1, false);
      c.w.a.b.resize(c.w.a.p.outputs());
      c.w.a.bnorm.resize(1);
      CompactOuter outer(c.w.a, workers);
      outer.known_solution = false;
      std::cout << std::setprecision(17);
      for (size_t q = 0; q < outer.x.size(); ++q) {
        auto z = chroma_bridge::source(q, chroma_bridge::probe_rhs);
        outer.x[q] = {z.re, z.im};
      }
      apply_cpu(c.w.a, outer.x.data(), outer.residual.data(), workers);
      std::cout << "STANDALONE original FP64 operator_match="
                << in.compare(outer.residual) << std::endl;
      std::string shader;
      for (auto name : {"kernels.metal", "bicgstab.metal", "clover.metal",
                        "multigrid.metal", "sparse_coarse.metal",
                        "compact_clover.metal", "mg_real_tune.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("missing shader");
        shader.append((std::istreambuf_iterator<char>(f)), {});
        shader += '\n';
      }
      MetalRun metal(shader);
      std::cout << "DEVICE " << [[metal.device name] UTF8String] << std::endl;
      std::ofstream trace(o.output + ".trace.csv"),
          summary(o.output + ".metal-solves.csv"),
          screens(o.output + ".screens.csv");
      trace << std::setprecision(17)
            << "method,rhs,outer,krylov,full_relative,fine_schur_applications,"
               "seconds\n";
      summary << std::setprecision(17)
              << "method,rhs,outer,krylov,seconds,residual_seconds,full_"
                 "relative,save_seconds,fine_schur_applications\n";
      screens << std::setprecision(17)
              << "method,pre,post,weight,bottom_limit,coarse_calls,coarse_"
                 "iterations,coarse_at_limit,bottom_iterations,bottom_max_"
                 "iterations,bottom_at_limit,bottom_worst_true,bottom_padded_"
                 "iterations,projected_seconds\n";
      CountedSchur op(metal, c);
      @autoreleasepool {
        check_schur(c, op);
      }
      release(c.links32);
      release(c.oddA32);
      release(c.inverseEven32);
      {
        @autoreleasepool {
          auto start = Clock::now();
          uint64_t setup_before = op.applications;
          Transfer first(c.parity.oe.p.rows, 12, o.vectors, h.shape, o.blocks,
                         true, workers, &c.parity.odd);
          build_gpu_filtered_basis(first, op, c, o.sweeps, o.filter_omega,
                                   workers);
          std::cout << "BASIS seconds=" << seconds(start)
                    << " filter_fine_S=" << op.applications - setup_before
                    << std::endl;
          BlockMatrix a1(2 * o.vectors, first.aggregates, workers);
          assemble_schur(a1, c, first);
          release(c.inverseEven);
          std::cout << "COARSE1 seconds=" << seconds(start)
                    << " entries=" << a1.col.size() << " dofs=" << a1.n
                    << std::endl;
          Transfer second(a1.rows, a1.bs, o.second_vectors, coarse_shape,
                          hypercube(2), false, workers);
          second.build(
              [&](const D *x, D *y) { a1.apply(x, y, o.second_vectors); },
              [&](const D *x, D *y) { a1.apply(x, y, o.second_vectors, true); },
              16);
          BlockMatrix a2(2 * o.second_vectors, second.aggregates, workers);
          project_matrix(a2, a1, second);
          first.release64();
          second.release64();
          a1.release64();
          a2.release64();
          std::cout << "HIERARCHY seconds=" << seconds(start)
                    << " coarse1=" << a1.n << " coarse2=" << a2.n << std::endl;
          TunedGPUCycle pre(first, metal, 1, c.inverseOdd32);
          release(first.basis32);
          release(c.inverseOdd32);
          pre.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                            id<MTLBuffer> y) { op.encode(e, x, y); };
          GPUBlockCoarseSolve<BlockMatrix> matrix1(a1, metal, 1, 1, .1),
              bottom(a2, metal, 1, o.bottom_limit, o.bottom_tol);
          GPUBlockOperator op1(matrix1);
          TunedGPUCycle pre1(second, metal, 1, a1.di32);
          release(second.basis32);
          pre1.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                             id<MTLBuffer> y) { op1.encode(e, x, y); };
          pre1.coarse = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                            id<MTLBuffer> y) { bottom.encode(e, x, y); };
          uint64_t bottom_at_limit = 0;
          int bottom_max = 0;
          pre1.check = [&]() {
            bottom.check();
            int it = int(static_cast<Cx<float> *>([bottom.q contents])[7].re);
            bottom_max = std::max(bottom_max, it);
            bottom_at_limit += it == bottom.limit;
          };
          GPUGenericFGMRES<GPUBlockOperator, TunedGPUCycle> coarse(a1.n, 1, op1,
                                                                   pre1, 8);
          coarse.allow_inexact = true;
          uint64_t coarse_calls = 0, coarse_iterations = 0, coarse_at_limit = 0;
          pre.blocking_coarse = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
            auto info = coarse.run_device(x, o.coarse_tol, o.coarse_limit);
            ++coarse_calls;
            coarse_iterations += info.max_iterations;
            coarse_at_limit += info.max_iterations == o.coarse_limit;
            pre.copy(coarse.delta, y, a1.n);
          };
          release(a1.a32);
          release(a1.di32);
          release(a2.a32);
          release(a2.di32);
          GPUGenericFGMRES<CountedSchur, TunedGPUCycle> mg(first.rows * 12, 1,
                                                           op, pre, o.restart);
          mg.allow_inexact = true;
          mg.blocking_pre = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
            pre.apply_device(x, y);
          };
          mg.progress = [](int k, double r) {
            if (k % 96 == 0)
              std::cout << "INNER metal-mg K=" << k << " relative_FP32=" << r
                        << std::endl;
          };
          CompactDeviceCorrection<decltype(mg)> correction(c, op, mg, workers);
          std::cout << "MG_READY setup_seconds=" << seconds(start)
                    << " allocated_GiB="
                    << double(metal.device.currentAllocatedSize) / (1ull << 30)
                    << std::endl;
          auto variants = cycles_arg(o.cycles, o.bottom_limit);
          size_t best = 0;
          double best_score = o.screen_outers ? 1e100 : 0;
          for (size_t i = 0; o.screen_outers && i < variants.size(); ++i) {
            auto v = variants[i];
            pre.pre_steps = v.pre;
            pre.post_steps = v.post;
            pre.omega = v.omega;
            bottom.limit = v.bottom;
            bottom.bp.limit = v.bottom;
            bottom.reset();
            bottom_at_limit = 0;
            bottom_max = 0;
            std::cout << "CYCLE method=" << v.name << " pre=" << v.pre
                      << " post=" << v.post << " omega=" << v.omega
                      << " bottom_cap=" << v.bottom << std::endl;
            source(c, outer, 0);
            uint64_t calls = coarse_calls, it = coarse_iterations,
                     cap = coarse_at_limit, bottom_it = bottom.iterations;
            auto r =
                refine(c, outer, correction, v.name, o.output, trace, summary,
                       0, op, o.screen_outers, o.screen_limit, false);
            double ratio = r.rho / r.previous,
                   extra =
                       (ratio > 0 && ratio < 1)
                           ? std::max(0., std::ceil(std::log(target / r.rho) /
                                                    std::log(ratio)))
                           : 1e6;
            double score = r.elapsed + extra * r.last_seconds;
            screens << v.name << ',' << v.pre << ',' << v.post << ',' << v.omega
                    << ',' << v.bottom << ',' << coarse_calls - calls << ','
                    << coarse_iterations - it << ',' << coarse_at_limit - cap
                    << ',' << bottom.iterations - bottom_it << ',' << bottom_max
                    << ',' << bottom_at_limit << ',' << bottom.worst_true << ','
                    << bottom.padded_iterations << ',' << score << '\n';
            screens.flush();
            std::cout << "PROJECTION method=" << v.name << " seconds=" << score
                      << " coarse_at_limit=" << coarse_at_limit - cap << "/"
                      << coarse_calls - calls << std::endl;
            if (score < best_score) {
              best_score = score;
              best = i;
            }
          }
          auto v = variants[best];
          std::cout << "SELECTED method=" << v.name << " pre=" << v.pre
                    << " post=" << v.post << " weight=" << v.omega
                    << " bottom_cap=" << v.bottom
                    << " projected_seconds=" << best_score << std::endl;
          if (o.finish) {
            pre.pre_steps = v.pre;
            pre.post_steps = v.post;
            pre.omega = v.omega;
            bottom.limit = v.bottom;
            bottom.bp.limit = v.bottom;
            for (int rhs = 0; rhs < o.rhs; ++rhs) {
              source(c, outer, rhs);
              refine(c, outer, correction, "metal-mg", o.output, trace, summary,
                     rhs, op, 12, 8192, true);
            }
          }
        }
      } // Release the hierarchy and MG workspaces before the BiCGStab timing.
      if (o.bicg) {
        std::cout << "BICG_READY allocated_GiB="
                  << double(metal.device.currentAllocatedSize) / (1ull << 30)
                  << std::endl;
        GPUBiCGStab bi(metal, c.parity.oe, 1, workers);
        bi.encode_override = [&](id<MTLComputeCommandEncoder> e,
                                 id<MTLBuffer> x,
                                 id<MTLBuffer> y) { op.encode(e, x, y); };
        CompactDeviceCorrection<decltype(bi)> bi_correction(c, op, bi, workers);
        for (int rhs = 0; rhs < o.rhs; ++rhs) {
          source(c, outer, rhs);
          refine(c, outer, bi_correction, "metal-bicgstab", o.output, trace,
                 summary, rhs, op, 12, 8192, true);
        }
      }
      rusage usage{};
      getrusage(RUSAGE_SELF, &usage);
      std::cout << "PEAK_RSS_GiB=" << double(usage.ru_maxrss) / (1ull << 30)
                << std::endl;
      std::cout << "PASS tuning complete; only SOLVED rows certify full FP64 "
                   "tolerance; native Chroma verification required"
                << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
    return 0;
  }
}
