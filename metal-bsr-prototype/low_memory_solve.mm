#include "low_memory_clover.hpp"
#include "mg_real_tune.hpp"
#include "near_null_svd.hpp"

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>

using namespace MetalCloverEnv;
namespace {
  constexpr int workers = 8;
  constexpr double tolerance = 1e-12;
  template <class T> void release(std::vector<T> &v) {
    std::vector<T>().swap(v);
  }
  struct CountedSchur : GPUCompactSchur {
    using GPUCompactSchur::GPUCompactSchur;
    uint64_t applications = 0;
    void encode(id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                id<MTLBuffer> y) {
      ++applications;
      GPUCompactSchur::encode(e, x, y);
    }
  };
  LatticeShape shape_arg(std::string v) {
    std::stringstream s(v);
    LatticeShape b{};
    for (int i = 0; i < 4; ++i) {
      char c;
      if (i && (!(s >> c) || c != ','))
        throw std::runtime_error("blocks syntax");
      if (!(s >> b[i]) || !b[i])
        throw std::runtime_error("block extent");
    }
    return b;
  }
  void fill_source(MappedField &f, const LowMemoryClover &c, uint32_t tag,
                   double scale) {
    auto *values = static_cast<D *>(f.data);
    uint32_t slab = c.w.stride(3) * 8;
    for (uint32_t b = 0; b < c.w.a.p.rows; b += slab) {
      uint32_t count = std::min(slab, c.w.a.p.rows - b);
      chunks(count * 12, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t i = size_t(b) * 12 + first; i < size_t(b) * 12 + last;
             ++i) {
          auto z = chroma_bridge::source(i, tag);
          values[i] = {scale * z.re, scale * z.im};
        }
      });
      f.drop(size_t(b) * 12 * 16, size_t(count) * 12 * 16, true);
    }
  }
  template <class Correction>
  CloverSolveInfo solve(LowMemoryClover &c, Correction &correction,
                        CountedSchur &op, MappedField &residual,
                        const std::string &output, const std::string &method,
                        int rhs, std::ostream &trace, std::ostream &csv) {
    uint32_t tag = 20261107 + rhs;
    size_t bytes = size_t(c.w.a.p.outputs()) * 16;
    std::string work = output + ".work-x.bin";
    MappedField x(work, bytes, true);
    auto *prepared = static_cast<Cx<float> *>([correction.input contents]);
    double scale = 1 / std::sqrt(double(c.w.a.p.outputs()));
    chunks(c.w.a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t i = first; i < last; ++i) {
        auto z = chroma_bridge::source(i, tag);
        prepared[i] = {float(scale * z.re), float(scale * z.im)};
      }
    });
    c.drop_coefficients();
    auto start = Clock::now();
    auto work_before = op.m.work;
    uint64_t before = op.applications;
    double rho = 1, check_s = 0;
    int k = 0, outer = 0;
    uint32_t slab = c.w.stride(3) * 8;
    for (outer = 1; outer <= 12; ++outer) {
      auto step = Clock::now();
      auto prev = op.applications;
      auto info = correction.run_prepared(.002, 16384);
      k += info.max_iterations;
      const auto *dx = correction.correction();
      auto *values = static_cast<D *>(x.data);
      residual.drop();
      for (uint32_t begin = 0; begin < c.w.a.p.rows; begin += slab) {
        uint32_t count = std::min(slab, c.w.a.p.rows - begin);
        chunks(count * 12, workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = size_t(begin) * 12 + first;
               i < size_t(begin) * 12 + last; ++i) {
            values[i].re += double(dx[i].re);
            values[i].im += double(dx[i].im);
          }
        });
        x.drop(size_t(begin) * 12 * 16, size_t(count) * 12 * 16, true);
      }
      auto check = Clock::now();
      rho = streamed_residual(c, x, residual, tag, workers, false, prepared);
      check_s += seconds(check);
      if (!std::isfinite(rho))
        throw std::runtime_error("nonfinite full FP64 residual");
      trace << method << ',' << rhs << ',' << outer << ','
            << info.max_iterations << ',' << rho << ','
            << op.applications - prev << ',' << seconds(step) << '\n';
      trace.flush();
      std::cout << "OUTER method=" << method << " rhs=" << rhs
                << " step=" << outer << " K=" << info.max_iterations
                << " fine_S=" << op.applications - prev
                << " seconds=" << seconds(step) << " full_FP64_relative=" << rho
                << std::endl;
      if (rho <= tolerance)
        break;
    }
    if (rho > tolerance)
      throw std::runtime_error("full FP64 tolerance failed");
    double elapsed = seconds(start);
    auto work_count = op.m.work - work_before;
    x.rename_to(output + "." + method + ".rhs" + std::to_string(rhs) + ".bin");
    csv << method << ',' << rhs << ',' << outer << ',' << k << ',' << elapsed
        << ',' << check_s << ',' << rho << ',' << op.applications - before
        << ',' << work_count.total() << ',' << work_count.gflops(elapsed) << ','
        << work_count.fine << ',' << work_count.coarse << ','
        << work_count.transfer << ',' << work_count.local_vector << '\n';
    csv.flush();
    std::cout << "SOLVED method=" << method << " rhs=" << rhs
              << " seconds=" << elapsed << " K=" << k
              << " fine_S=" << op.applications - before
              << " full_FP64_relative=" << rho
              << " nominal_operator_GFLOP_per_s=" << work_count.gflops(elapsed)
              << std::endl;
    return {outer, k,         op.applications - before, elapsed, check_s,
            rho,   work_count};
  }
} // namespace
int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      if (argc == 2 && std::string(argv[1]) == "--basis-self-test") {
        check_near_null_svd();
        return 0;
      }
      if (argc < 8)
        throw std::runtime_error(
            "Usage: low_memory_solve BRIDGE_PREFIX OUTPUT_PREFIX mg|bicg RHS "
            "BLOCKS RESTART FILTER_SWEEPS [SECOND_FILTER_SWEEPS BOTTOM_LIMIT] "
            "[--basis filtered|near-null-svd --null-vectors N --null-iters N "
            "--null-restart N --null-warm N --null-tol X --null-batch N]");
      std::string prefix = argv[1], output = argv[2], method = argv[3];
      int rhs = std::stoi(argv[4]), restart = std::stoi(argv[6]),
          sweeps = std::stoi(argv[7]);
      auto blocks = shape_arg(argv[5]);
      int second_sweeps = 16, bottom_limit = 4, first_option = 8;
      if (argc >= 10 && std::string(argv[8]).rfind("--", 0) != 0) {
        second_sweeps = std::stoi(argv[8]);
        bottom_limit = std::stoi(argv[9]);
        first_option = 10;
      }
      std::string basis = "filtered";
      NearNullOptions null_options;
      for (int i = first_option; i < argc; ++i) {
        std::string key = argv[i];
        if (++i >= argc)
          throw std::runtime_error("missing basis option value");
        std::string value = argv[i];
        if (key == "--basis")
          basis = value;
        else if (key == "--null-vectors")
          null_options.candidates = std::stoi(value);
        else if (key == "--null-iters")
          null_options.iterations = std::stoi(value);
        else if (key == "--null-restart")
          null_options.restart = std::stoi(value);
        else if (key == "--null-warm")
          null_options.warm = std::stoi(value);
        else if (key == "--null-batch")
          null_options.batch = std::stoi(value);
        else if (key == "--null-tol")
          null_options.tolerance = std::stod(value);
        else
          throw std::runtime_error("unknown basis option " + key);
      }
      if (basis != "filtered" && basis != "near-null-svd")
        throw std::runtime_error("unknown basis construction");
      if ((method != "mg" && method != "bicg") || rhs < 1 || rhs > 192 ||
          restart < 2 || restart > 12 || sweeps < 1 || second_sweeps < 1 ||
          bottom_limit < 1)
        throw std::runtime_error("option range");
      if (std::filesystem::exists(output + "." + method + ".csv"))
        throw std::runtime_error("output exists");
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
      std::ifstream input(prefix + ".bridge", std::ios::binary);
      auto h = production::read_header(input);
      size_t offset = size_t(input.tellg());
      input.close();
      MappedField bridge(prefix + ".bridge");
      if (bridge.bytes != offset + size_t(h.volume()) * 48 * 16)
        throw std::runtime_error("bridge size");
      LowMemoryClover c(bridge, h, offset, workers);
      std::cout << std::setprecision(17);
      size_t bytes = size_t(c.w.a.p.outputs()) * 16;
      {
        MappedField x(output + ".probe-x.bin", bytes, true), r("", bytes, true);
        fill_source(x, c, chroma_bridge::probe_rhs, 1);
        double error = streamed_residual(c, x, r, 0, workers, true);
        std::cout << "FP64 original_Chroma_operator_match=" << error
                  << std::endl;
        if (!std::isfinite(error) || error > 1e-13)
          throw std::runtime_error("original Chroma full operator mismatch");
      }
      std::filesystem::remove(output + ".probe-x.bin");
      std::string shader;
      for (auto name :
           {"kernels.metal", "bicgstab.metal", "clover.metal",
            "multigrid.metal", "sparse_coarse.metal", "compact_clover.metal",
            "mg_real_tune.metal", "packed_clover.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("shader missing");
        shader.append((std::istreambuf_iterator<char>(f)), {});
        shader += '\n';
      }
      MetalRun metal(shader);
      CountedSchur op(metal, c);
      release(c.links32);
      release(c.oddA32);
      release(c.inverseEven32);
      c.drop_coefficients();
      MappedField residual("", bytes, true);
      std::ofstream trace(output + "." + method + ".trace.csv"),
          csv(output + "." + method + ".csv");
      trace << std::setprecision(17)
            << "method,rhs,outer,krylov,full_relative,fine_schur_applications,"
               "seconds\n";
      csv << std::setprecision(17)
          << "method,rhs,outer,krylov,seconds,residual_seconds,full_relative,"
             "fine_schur_applications,nominal_operator_flops,operator_gflops,"
             "fine_flops,coarse_flops,transfer_flops,local_vector_flops\n";
      std::vector<CloverSolveInfo> solved;
      if (method == "bicg") {
        GPUBiCGStab bi(metal, c.parity.oe, 1, workers);
        bi.encode_override = [&](id<MTLComputeCommandEncoder> e,
                                 id<MTLBuffer> x,
                                 id<MTLBuffer> y) { op.encode(e, x, y); };
        bi.progress = [](int k, double r) {
          if (k % 256 == 0)
            std::cout << "INNER BiCGStab K=" << k << " relative=" << r
                      << std::endl;
        };
        CompactDeviceCorrection<decltype(bi)> correction(c, op, bi, workers);
        std::cout << "BICG_READY GPU_GiB="
                  << double(metal.device.currentAllocatedSize) / (1ull << 30)
                  << std::endl;
        for (int i = 0; i < rhs; ++i)
          solved.push_back(solve(c, correction, op, residual, output,
                                 "metal-bicgstab", i, trace, csv));
      } else {
        auto start = Clock::now();
        auto before = op.applications;
        LatticeShape coarse_shape{};
        for (int mu = 0; mu < 4; ++mu) {
          if (h.shape[mu] % blocks[mu] || (h.shape[mu] / blocks[mu]) % 2)
            throw std::runtime_error("aggregation geometry");
          coarse_shape[mu] = h.shape[mu] / blocks[mu];
        }
        Transfer first(c.parity.oe.p.rows, 12, 8, h.shape, blocks, true,
                       workers, &c.parity.odd);
        if (basis == "near-null-svd") {
          std::ofstream null_log(output + ".near-null.csv"),
              spectrum(output + ".svd.csv");
          null_log << std::setprecision(17);
          spectrum << std::setprecision(17);
          build_gpu_near_null_svd(first, op, null_options, null_log, spectrum);
        } else {
          build_gpu_filtered_basis(first, op, c, sweeps, .65, workers);
        }
        c.drop_coefficients();
        std::cout << "BASIS method=" << basis << " seconds=" << seconds(start)
                  << " fine_S=" << op.applications - before << std::endl;
        BlockMatrix a1(16, first.aggregates, workers);
        assemble_schur(a1, c, first);
        c.drop_coefficients();
        std::cout << "COARSE1 seconds=" << seconds(start) << " dofs=" << a1.n
                  << std::endl;
        Transfer second(a1.rows, a1.bs, 8, coarse_shape, hypercube(2), false,
                        workers);
        second.build([&](const D *x, D *y) { a1.apply(x, y, 8); },
                     [&](const D *x, D *y) { a1.apply(x, y, 8, true); },
                     second_sweeps);
        BlockMatrix a2(16, second.aggregates, workers);
        project_matrix(a2, a1, second);
        first.release64();
        second.release64();
        a1.release64();
        a2.release64();
        std::vector<Cx<float>> empty;
        TunedGPUCycle pre(first, metal, 1, empty, 2);
        pre.inv = op.Aodd;
        release(first.basis32);
        pre.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                          id<MTLBuffer> y) { op.encode(e, x, y); };
        GPUBlockCoarseSolve<BlockMatrix> matrix1(a1, metal, 1, 1, .1),
            bottom(a2, metal, 1, bottom_limit, .1);
        GPUBlockOperator op1(matrix1);
        TunedGPUCycle pre1(second, metal, 1, a1.di32);
        release(second.basis32);
        pre1.product = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                           id<MTLBuffer> y) { op1.encode(e, x, y); };
        pre1.coarse = [&](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                          id<MTLBuffer> y) { bottom.encode(e, x, y); };
        pre1.check = [&]() { bottom.check(); };
        GPUGenericFGMRES<GPUBlockOperator, TunedGPUCycle> coarse(a1.n, 1, op1,
                                                                 pre1, 8);
        coarse.allow_inexact = true;
        uint64_t coarse_calls = 0, coarse_iterations = 0, coarse_at_limit = 0;
        pre.blocking_coarse = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
          auto info = coarse.run_device(x, .1, 32);
          ++coarse_calls;
          coarse_iterations += info.max_iterations;
          coarse_at_limit += info.max_iterations == 32;
          pre.copy(coarse.delta, y, a1.n);
        };
        release(a1.a32);
        release(a1.di32);
        release(a2.a32);
        release(a2.di32);
        GPUGenericFGMRES<CountedSchur, TunedGPUCycle> mg(first.rows * 12, 1, op,
                                                         pre, restart);
        mg.blocking_pre = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
          pre.apply_device(x, y);
        };
        mg.progress = [](int k, double r) {
          std::cout << "INNER MG K=" << k << " relative=" << r << std::endl;
        };
        CompactDeviceCorrection<decltype(mg)> correction(c, op, mg, workers);
        std::cout << "MG_READY hierarchy_seconds=" << seconds(start)
                  << " hierarchy_fine_S=" << op.applications - before
                  << " GPU_GiB="
                  << double(metal.device.currentAllocatedSize) / (1ull << 30)
                  << " recommended_GiB="
                  << double(metal.device.recommendedMaxWorkingSetSize) /
                         (1ull << 30)
                  << std::endl;
        for (int i = 0; i < rhs; ++i)
          solved.push_back(solve(c, correction, op, residual, output,
                                 "metal-mg", i, trace, csv));
        std::cout << "COARSE_STATS calls=" << coarse_calls
                  << " iterations=" << coarse_iterations
                  << " at_limit=" << coarse_at_limit
                  << " bottom_iterations=" << bottom.iterations
                  << " bottom_padded=" << bottom.padded_iterations
                  << " bottom_worst_true=" << bottom.worst_true << std::endl;
      }
      double solve_seconds = 0;
      uint64_t solve_flops = 0;
      for (const auto &info : solved) {
        solve_seconds += info.seconds;
        solve_flops += info.operator_work.total();
      }
      std::cout << "TOTAL_SOLVE method=" << method << " rhs=" << solved.size()
                << " seconds=" << solve_seconds
                << " mean_seconds=" << solve_seconds / solved.size()
                << " nominal_operator_flops=" << solve_flops
                << " nominal_operator_GFLOP_per_s="
                << double(solve_flops) / (solve_seconds * 1e9) << std::endl;
      rusage u{};
      getrusage(RUSAGE_SELF, &u);
      std::cout << "PASS all full FP64 residuals <=1e-12; original Chroma "
                   "final verification required; peak_RSS_GiB="
                << double(u.ru_maxrss) / (1ull << 30) << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "ERROR " << e.what() << std::endl;
      return 1;
    }
    return 0;
  }
}
