#include "sparse_checks.hpp"
#include <memory>
#include <sys/resource.h>

#include "sparse_rhs.hpp"

using namespace MetalCloverEnv;

template <class Inner, class Pre, class Coarse, class FG>
void benchmark_sparse(CloverSystem &c, const TwoLevel &h,
                      const SparseCoarse &ac, Inner &correction, Pre &pre,
                      Coarse &coarse, FG &fg, int workers, uint32_t total_rhs,
                      int samples, bool warmup, bool self,
                      const std::string &name, double workspace_ms,
                      double problem_ms, double tol, double inner_tol,
                      std::ofstream &csv, std::ofstream &trace) {
  uint32_t tile = c.w.a.p.rhs;
  OuterWorkspace outer(c.w.a, workers);
  auto once = [&](uint32_t offset) {
    set_rhs(c, offset, workers, self);
    coarse.reset();
    pre.calls = pre.schur_calls = fg.fine_calls = 0;
    correction.reset_metrics();
    return solve(c.w.a, correction, outer, workers, tol, inner_tol, 256);
  };
  if (warmup && !self)
    once(0);
  for (int trial = 0; trial < samples; ++trial) {
    double total = 0, inner = 0, residual = 0, update = 0, prep = 0, recon = 0,
           worst = 0, error = 0, worst_coarse = 0, rhs_setup = 0;
    uint64_t fine = 0, pre_calls = 0, coarse_ops = 0, coarse_iters = 0;
    int outer_min = 100, outer_max = 0, krylov_sum = 0, krylov_min = 10000,
        krylov_max = 0;
    uint32_t tiles = 0;
    for (uint32_t offset = 0; offset < total_rhs; offset += tile) {
      auto start = Clock::now();
      set_rhs(c, offset, workers, self);
      rhs_setup += seconds(start) * 1000;
      coarse.reset();
      pre.calls = pre.schur_calls = fg.fine_calls = 0;
      correction.reset_metrics();
      auto r = solve(c.w.a, correction, outer, workers, tol, inner_tol, 256);
      ++tiles;
      total += r.total_ms;
      inner += r.inner_ms;
      residual += r.residual_ms;
      update += r.update_ms;
      prep += correction.prepare_ms;
      recon += correction.reconstruct_ms;
      worst = std::max(worst, r.max_residual);
      error = std::max(error, r.solution_error);
      fine += pre.schur_calls + fg.fine_calls;
      pre_calls += pre.calls;
      coarse_ops += coarse.applications;
      coarse_iters += coarse.iterations;
      worst_coarse = std::max(worst_coarse, coarse.worst_true);
      outer_min = std::min(outer_min, r.outer);
      outer_max = std::max(outer_max, r.outer);
      krylov_sum += r.inner_max_sum;
      krylov_min = std::min(krylov_min, r.inner_max_sum);
      krylov_max = std::max(krylov_max, r.inner_max_sum);
      std::cout << "L=" << h.extent << " mass=" << c.w.a.mass
                << " RHS=" << total_rhs << " tile=" << offset << "+" << tile
                << " " << name << " trial=" << trial << " ms=" << r.total_ms
                << " outer=" << r.outer << " inner=" << r.inner_max_sum
                << " fineS=" << pre.schur_calls + fg.fine_calls
                << " coarse_spmv=" << coarse.applications
                << " coarse_true=" << coarse.worst_true
                << " true_r=" << r.max_residual << " PASS" << std::endl;
      if (trace)
        for (auto tr : r.trace)
          trace << std::setprecision(12) << h.extent << ',' << c.w.a.mass << ','
                << total_rhs << ',' << tile << ',' << offset << ',' << name
                << ',' << trial << ',' << tr.outer << ',' << tr.max_residual
                << ',' << tr.inner.min_iterations << ','
                << tr.inner.max_iterations << '\n';
    }
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    if (csv) {
      csv << std::setprecision(12) << h.extent << ',' << total_rhs << ','
          << tile << ',' << tiles << ',' << c.w.a.mass << ',' << c.coefficient
          << ',' << h.nvec << ',' << ac.n << ',' << ac.bs << ','
          << ac.col.size() << ',' << h.sweeps << ',' << h.omega << ','
          << fg.restart << ',' << coarse.limit << ',' << coarse.tol << ','
          << name << ',' << trial << ',' << tol << ',' << total << ',' << inner
          << ',' << residual << ',' << update << ',' << prep << ',' << recon
          << ',' << h.setup_ms + ac.setup_ms << ',' << workspace_ms << ','
          << problem_ms << ',' << rhs_setup << ',' << outer_min << ','
          << outer_max << ',' << krylov_sum << ',' << krylov_min << ','
          << krylov_max << ',' << fine << ',' << pre_calls << ',' << coarse_ops
          << ',' << coarse_iters << ',' << worst_coarse << ',' << worst << ','
          << error << ',' << total / total_rhs << ','
          << 1000 * total_rhs / total << ','
          << double(usage.ru_maxrss) / 1048576 << '\n';
      csv.flush();
      trace.flush();
    }
    std::cout << "BATCH " << name << " L=" << h.extent << " mass=" << c.w.a.mass
              << " RHS=" << total_rhs << " tile=" << tile << " trial=" << trial
              << " total_ms=" << total << " true_r=" << worst
              << " peak_RSS_MiB=" << double(usage.ru_maxrss) / 1048576
              << std::endl;
  }
}

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      std::vector<uint32_t> sizes{16}, rhs_values{96, 192};
      std::vector<double> masses{.1, .5};
      uint32_t tile = 96, nvec = 12, block = 4, restart = 8;
      int workers = 8, sweeps = 8, samples = 3, coarse_limit = 4;
      double omega = 1, coarse_tol = .1, tol = 1e-10, inner_tol = .002;
      bool self = false, validate = false, warmup = true, smooth_only = false;
      std::string backend = "both", csv_path;
      for (int j = 1; j < argc; ++j) {
        std::string opt = argv[j];
        auto next = [&]() {
          if (++j >= argc)
            throw std::runtime_error("missing value");
          return std::string(argv[j]);
        };
        if (opt == "--sizes")
          sizes = list_arg(next(), 16);
        else if (opt == "--rhs")
          rhs_values = list_arg(next(), 192);
        else if (opt == "--tile-rhs")
          tile = uint32_t(std::stoul(next()));
        else if (opt == "--vectors")
          nvec = uint32_t(std::stoul(next()));
        else if (opt == "--aggregate")
          block = uint32_t(std::stoul(next()));
        else if (opt == "--restart")
          restart = uint32_t(std::stoul(next()));
        else if (opt == "--threads")
          workers = std::stoi(next());
        else if (opt == "--setup-sweeps")
          sweeps = std::stoi(next());
        else if (opt == "--samples")
          samples = std::stoi(next());
        else if (opt == "--omega")
          omega = std::stod(next());
        else if (opt == "--coarse-iters")
          coarse_limit = std::stoi(next());
        else if (opt == "--coarse-tol")
          coarse_tol = std::stod(next());
        else if (opt == "--backend")
          backend = next();
        else if (opt == "--csv")
          csv_path = next();
        else if (opt == "--mass") {
          masses.clear();
          std::stringstream in(next());
          std::string s;
          while (std::getline(in, s, ','))
            masses.push_back(std::stod(s));
        } else if (opt == "--self-test")
          self = true;
        else if (opt == "--validate-only")
          validate = true;
        else if (opt == "--no-warmup")
          warmup = false;
        else if (opt == "--smoother-only")
          smooth_only = true;
        else if (opt == "--csw") {
          if (std::stod(next()) != 1)
            throw std::runtime_error("this comparison keeps c_sw=1");
        } else if (opt == "--help") {
          std::cout
              << "--sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --vectors "
                 "12 --aggregate 4 --restart 8 --coarse-iters 4 --coarse-tol "
                 ".1 --samples 3 --threads 8 --csw 1 --backend both|cpu|metal "
                 "--csv results/sparse-mg.csv --self-test --validate-only "
                 "--no-warmup --smoother-only\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + opt);
      }
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
      if (self) {
        sizes = {2, 4};
        rhs_values = {3};
        tile = 3;
        nvec = 6;
        block = 2;
        samples = 1;
        restart = 6;
      }
      if (tile < 1 || tile > 96 || nvec < 1 || nvec > 24 || block < 2 ||
          block % 2 || workers < 1 || workers > 32 || sweeps < 0 ||
          sweeps > 96 || restart < 2 || restart > 12 || samples < 1 ||
          samples > 9 || samples % 2 == 0 || coarse_limit < 1 ||
          coarse_limit > 64 || !(coarse_tol > 0 && coarse_tol < 1) ||
          !(omega > 0 && omega <= 1.5) || masses.empty() ||
          (backend != "both" && backend != "cpu" && backend != "metal"))
        throw std::runtime_error("invalid options");
      for (auto m : masses)
        if (!(m > 0 && m <= 2))
          throw std::runtime_error("invalid mass");
      for (auto l : sizes)
        if (l % block || l % 2)
          throw std::runtime_error("aggregate must divide even extent");
      for (auto rhs : rhs_values)
        if (rhs % tile)
          throw std::runtime_error("tile-rhs must divide total RHS");
      std::string source;
      for (auto name : {"kernels.metal", "bicgstab.metal", "clover.metal",
                        "multigrid.metal", "sparse_coarse.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("run from prototype directory");
        source.append((std::istreambuf_iterator<char>(f)), {});
        source += '\n';
      }
      MetalRun metal(source);
      std::cout << "GPU=" << [[metal.device name] UTF8String]
                << " maxBuffer_MiB=" << metal.device.maxBufferLength / 1048576
                << " sparse coarse BiCGStab; csw=1; physical tile=" << tile
                << " workers=" << workers << std::endl;
      std::ofstream csv, trace;
      if (!csv_path.empty()) {
        csv.open(csv_path);
        trace.open(csv_path + ".trace.csv");
        if (!csv || !trace)
          throw std::runtime_error("cannot create results");
        csv << "extent,rhs,tile_rhs,tiles,mass,csw,vectors,coarse_dim,coarse_"
               "block,coarse_nnz_blocks,setup_sweeps,omega,restart,coarse_"
               "limit,coarse_tol,backend,trial,tolerance,total_ms,inner_ms,"
               "residual_ms,update_ms,prepare_ms,reconstruct_ms,hierarchy_"
               "setup_ms,workspace_setup_ms,problem_setup_ms,rhs_setup_ms,"
               "outer_min,outer_max,krylov_sum,krylov_min,krylov_max,fine_"
               "schur_calls,preconditioner_calls,coarse_spmv_calls,coarse_"
               "iterations,worst_coarse_true_residual,max_true_residual,"
               "solution_error,ms_per_rhs,rhs_per_second,peak_rss_mib\n";
        trace << "extent,mass,rhs,tile_rhs,offset,backend,trial,outer_step,max_"
                 "true_residual,inner_min_iterations,inner_max_iterations\n";
      }
      int case_count = 0;
      for (auto mass : masses)
        for (auto l : sizes) {
          int first_backend = case_count++ % 2;
          TwoLevel h(l, {block, block, block, block}, nvec, mass, .25, 1, true,
                     workers, sweeps, omega, .65, 1, 1, false);
          SparseCoarse ac(h, mass, .25, 1, true, workers);
          auto t = Clock::now();
          CloverSystem c(l, tile, mass, .25, 1, true, workers, self);
          std::vector<Z>().swap(c.w.a.p.x);
          std::vector<Z>().swap(c.parity.eo.p.x);
          std::vector<Z>().swap(c.parity.oe.p.x);
          double problem_ms = seconds(t) * 1000;
          if (self || validate)
            check_sparse(h, ac, c, metal, workers, h.coarse <= 512);
          if (validate)
            continue;
          for (int slot = 0; slot < 2; ++slot) {
            int kind = (slot + first_backend) % 2;
            if ((kind == 0 && backend == "metal") ||
                (kind == 1 && backend == "cpu"))
              continue;
            @autoreleasepool {
              t = Clock::now();
              if (kind == 0) {
                CPUCloverSchur<float> op(c, workers);
                CPUMultigrid mg(c, h, op, workers);
                CPUCoarseSolve coarse(ac, tile, workers, coarse_limit,
                                      coarse_tol);
                mg.coarse_enabled = !smooth_only;
                mg.coarse_override = [&](const Cx<float> *x, Cx<float> *y) {
                  coarse.solve(x, y);
                };
                CPUFGMRES fg(c, op, mg, workers, restart);
                HostCloverCorrection<CPUFGMRES, float> correction(c, fg,
                                                                  workers);
                double setup = seconds(t) * 1000;
                for (auto rhs : rhs_values)
                  benchmark_sparse(
                      c, h, ac, correction, mg, coarse, fg, workers, rhs,
                      samples, warmup && rhs == rhs_values.front(), self,
                      smooth_only ? "cpu_smooth" : "cpu_sparse_mg", setup,
                      problem_ms, tol, inner_tol, csv, trace);
              } else {
                GPUCloverSchur op(metal, c);
                GPUMultigrid mg(c, h, op);
                GPUCoarseSolve coarse(ac, metal, tile, coarse_limit,
                                      coarse_tol);
                mg.coarse_enabled = !smooth_only;
                mg.coarse_override = [&](id<MTLComputeCommandEncoder> enc,
                                         id<MTLBuffer> x, id<MTLBuffer> y) {
                  coarse.encode(enc, x, y);
                };
                mg.coarse_check = [&]() {
                  if (mg.coarse_enabled)
                    coarse.check();
                };
                GPUFGMRES fg(c, op, mg, restart);
                DeviceCloverCorrectionT<GPUFGMRES> correction(c, op, fg,
                                                              workers);
                double setup = seconds(t) * 1000;
                std::cout << "Metal allocated_MiB="
                          << metal.device.currentAllocatedSize / 1048576
                          << std::endl;
                for (auto rhs : rhs_values)
                  benchmark_sparse(
                      c, h, ac, correction, mg, coarse, fg, workers, rhs,
                      samples, warmup && rhs == rhs_values.front(), self,
                      smooth_only ? "metal_smooth" : "metal_sparse_mg", setup,
                      problem_ms, tol, inner_tol, csv, trace);
              }
            }
          }
        }
      if (validate)
        std::cout << "PASS: sparse hierarchy and coarse-solver validation; no "
                     "full solves requested."
                  << std::endl;
      else
        std::cout << "PASS: all requested systems; full CPU FP64 clover-Wilson "
                     "residual checks for every RHS."
                  << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
