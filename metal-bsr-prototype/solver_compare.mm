#include "sparse_checks.hpp"
#include "sparse_rhs.hpp"

#include <memory>
#include <sys/resource.h>

using namespace MetalCloverEnv;

struct Work {
  uint64_t fine = 0, pre = 0, coarse_ops = 0, coarse_iters = 0;
  double coarse_true = 0;
};
struct Outcome {
  SolveResult result;
  Work work;
  double prepare = 0, reconstruct = 0;
};
struct Batch {
  double total = 0, inner = 0, residual = 0, update = 0, prepare = 0,
         reconstruct = 0, rhs_setup = 0, worst = 0, error = 0, coarse_true = 0;
  uint64_t fine = 0, pre = 0, coarse_ops = 0, coarse_iters = 0;
  int outer_min = 100, outer_max = 0, krylov_min = 10000, krylov_max = 0,
      krylov_sum = 0, dispatched = 0, tiles = 0;
  void add(const Outcome &s, double setup) {
    const auto &r = s.result;
    ++tiles;
    total += r.total_ms;
    inner += r.inner_ms;
    residual += r.residual_ms;
    update += r.update_ms;
    prepare += s.prepare;
    reconstruct += s.reconstruct;
    rhs_setup += setup;
    worst = std::max(worst, r.max_residual);
    error = std::max(error, r.solution_error);
    coarse_true = std::max(coarse_true, s.work.coarse_true);
    fine += s.work.fine;
    pre += s.work.pre;
    coarse_ops += s.work.coarse_ops;
    coarse_iters += s.work.coarse_iters;
    outer_min = std::min(outer_min, r.outer);
    outer_max = std::max(outer_max, r.outer);
    krylov_min = std::min(krylov_min, r.inner_max_sum);
    krylov_max = std::max(krylov_max, r.inner_max_sum);
    krylov_sum += r.inner_max_sum;
    dispatched += r.dispatched;
  }
};
struct Options {
  std::vector<uint32_t> sizes{16}, rhs{96, 192};
  std::vector<double> masses{.1, .5};
  uint32_t tile = 96, nvec = 12, block = 4, restart = 8;
  int workers = 8, samples = 3, bi_chunk = 1, coarse_limit = 4;
  double coarse_tol = .1;
  bool self = false, warmup = true, validate = false, metal_first = false;
  std::string csv, backend = "both", methods = "bicgstab,smooth,sparse_mg";
};

void compare_measure(CloverSystem &c, const TwoLevel &h, const SparseCoarse &ac,
                     OuterWorkspace &outer, const Options &o,
                     const std::string &method, const std::string &backend,
                     double inner_setup, double outer_setup,
                     double problem_setup, std::function<Outcome()> run,
                     std::ofstream &csv, std::ofstream &trace) {
  (void)outer;
  uint32_t tile = c.w.a.p.rhs;
  bool mg = method == "sparse_mg";
  std::string name = backend + "_" + method;
  if (o.warmup && !o.self) {
    set_rhs(c, 0, o.workers);
    run();
  }
  for (int trial = 0; trial < o.samples; ++trial) {
    Batch sum;
    for (uint32_t offset = 0; offset < o.rhs.back(); offset += tile) {
      auto t = Clock::now();
      set_rhs(c, offset, o.workers, o.self);
      double rhs_setup = seconds(t) * 1000;
      auto sample = run();
      sum.add(sample, rhs_setup);
      std::cout << "TILE L=" << h.extent << " mass=" << c.w.a.mass << " "
                << name << " trial=" << trial << " offset=" << offset
                << " ms=" << sample.result.total_ms
                << " outer=" << sample.result.outer
                << " krylov=" << sample.result.inner_max_sum
                << " fineS=" << sample.work.fine
                << " true_r=" << sample.result.max_residual << " PASS"
                << std::endl;
      if (trace)
        for (auto tr : sample.result.trace)
          trace << std::setprecision(12) << h.extent << ',' << c.w.a.mass << ','
                << tile << ',' << offset << ',' << name << ',' << trial << ','
                << tr.outer << ',' << tr.max_residual << ','
                << tr.inner.min_iterations << ',' << tr.inner.max_iterations
                << ',' << tr.inner.dispatched << '\n';
      uint32_t total_rhs = offset + tile;
      if (std::find(o.rhs.begin(), o.rhs.end(), total_rhs) == o.rhs.end())
        continue;
      rusage usage{};
      getrusage(RUSAGE_SELF, &usage);
      if (csv) {
        csv << std::setprecision(12) << h.extent << ',' << total_rhs << ','
            << tile << ',' << sum.tiles << ',' << c.w.a.mass << ",1," << name
            << ',' << trial << ",1e-10,.002," << o.workers << ',' << o.restart
            << ',' << o.bi_chunk << ',' << h.nvec << ',' << h.coarse << ','
            << ac.col.size() << ',' << o.coarse_limit << ',' << o.coarse_tol
            << ',' << sum.total << ',' << sum.inner << ',' << sum.residual
            << ',' << sum.update << ',' << sum.prepare << ',' << sum.reconstruct
            << ',' << (mg ? h.setup_ms + ac.setup_ms : 0) << ',' << inner_setup
            << ',' << outer_setup << ',' << problem_setup << ','
            << sum.rhs_setup << ',' << sum.outer_min << ',' << sum.outer_max
            << ',' << sum.krylov_min << ',' << sum.krylov_max << ','
            << sum.krylov_sum << ',' << sum.dispatched << ',' << sum.fine << ','
            << sum.pre << ',' << sum.coarse_ops << ',' << sum.coarse_iters
            << ',' << sum.coarse_true << ',' << sum.worst << ',' << sum.error
            << ',' << 1000 * total_rhs / sum.total << ','
            << double(usage.ru_maxrss) / 1048576 << '\n';
        csv.flush();
        trace.flush();
      }
      std::cout << "BATCH L=" << h.extent << " mass=" << c.w.a.mass
                << " RHS=" << total_rhs << " " << name << " trial=" << trial
                << " total_ms=" << sum.total << " true_r=" << sum.worst
                << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      Options o;
      for (int j = 1; j < argc; ++j) {
        std::string flag = argv[j];
        auto next = [&]() {
          if (++j >= argc)
            throw std::runtime_error("missing value");
          return std::string(argv[j]);
        };
        if (flag == "--sizes")
          o.sizes = list_arg(next(), 16);
        else if (flag == "--rhs")
          o.rhs = list_arg(next(), 192);
        else if (flag == "--tile-rhs")
          o.tile = uint32_t(std::stoul(next()));
        else if (flag == "--samples")
          o.samples = std::stoi(next());
        else if (flag == "--threads")
          o.workers = std::stoi(next());
        else if (flag == "--bi-chunk")
          o.bi_chunk = std::stoi(next());
        else if (flag == "--coarse-iters")
          o.coarse_limit = std::stoi(next());
        else if (flag == "--coarse-tol")
          o.coarse_tol = std::stod(next());
        else if (flag == "--backend")
          o.backend = next();
        else if (flag == "--methods")
          o.methods = next();
        else if (flag == "--csv")
          o.csv = next();
        else if (flag == "--mass") {
          o.masses.clear();
          std::stringstream in(next());
          std::string value;
          while (std::getline(in, value, ','))
            o.masses.push_back(std::stod(value));
        } else if (flag == "--csw") {
          if (std::stod(next()) != 1)
            throw std::runtime_error("comparison requires c_sw=1");
        } else if (flag == "--self-test")
          o.self = true;
        else if (flag == "--validate-only")
          o.validate = true;
        else if (flag == "--metal-first")
          o.metal_first = true;
        else if (flag == "--no-warmup")
          o.warmup = false;
        else if (flag == "--help") {
          std::cout
              << "--sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --samples "
                 "3 --threads 8 --csw 1 --bi-chunk 1 --coarse-iters 4 "
                 "--coarse-tol .1 --backend both|cpu|metal --methods "
                 "bicgstab,smooth,sparse_mg --csv results/solver-compare16.csv "
                 "--self-test --validate-only --metal-first --no-warmup\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + flag);
      }
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
      if (o.self) {
        o.sizes = {2, 4};
        o.rhs = {3, 6};
        o.tile = 3;
        o.nvec = 6;
        o.block = 2;
        o.restart = 6;
        o.samples = 1;
      }
      if (o.tile < 1 || o.tile > 96 || o.samples < 1 || o.samples > 9 ||
          o.samples % 2 == 0 || o.workers < 1 || o.workers > 32 ||
          o.bi_chunk < 1 || o.bi_chunk > 8 || o.coarse_limit < 1 ||
          o.coarse_limit > 64 || !(o.coarse_tol > 0 && o.coarse_tol < 1) ||
          (o.backend != "both" && o.backend != "cpu" && o.backend != "metal") ||
          o.masses.empty())
        throw std::runtime_error("invalid options");
      for (auto m : o.masses)
        if (!(m > -1 && m <= 2))
          throw std::runtime_error("bare mass must lie in (-1,2]");
      for (auto l : o.sizes)
        if (l % o.block || l % 2)
          throw std::runtime_error("aggregate must divide even extent");
      for (auto r : o.rhs)
        if (r % o.tile)
          throw std::runtime_error("tile width must divide total RHS");
      std::sort(o.rhs.begin(), o.rhs.end());
      if (std::adjacent_find(o.rhs.begin(), o.rhs.end()) != o.rhs.end())
        throw std::runtime_error("duplicate RHS");
      std::vector<std::string> methods;
      {
        std::stringstream in(o.methods);
        std::string s;
        while (std::getline(in, s, ',')) {
          if (s != "bicgstab" && s != "smooth" && s != "sparse_mg")
            throw std::runtime_error("invalid method");
          if (std::find(methods.begin(), methods.end(), s) != methods.end())
            throw std::runtime_error("duplicate method");
          methods.push_back(s);
        }
      }
      if (methods.empty())
        throw std::runtime_error("empty method list");
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
                << " CPU workers=" << o.workers
                << " c_sw=1 resident_RHS=" << o.tile
                << " BiCG_chunk=" << o.bi_chunk
                << " paired prefix batch timings" << std::endl;
      std::ofstream csv, trace;
      if (!o.csv.empty()) {
        csv.open(o.csv);
        trace.open(o.csv + ".trace.csv");
        if (!csv || !trace)
          throw std::runtime_error("cannot create results");
        csv << "extent,rhs,tile_rhs,tiles,mass,csw,backend,trial,tolerance,"
               "inner_tolerance,workers,restart,bi_chunk,vectors,coarse_dim,"
               "coarse_nnz_blocks,coarse_limit,coarse_tolerance,total_ms,inner_"
               "ms,residual_ms,update_ms,prepare_ms,reconstruct_ms,hierarchy_"
               "setup_ms,inner_workspace_setup_ms,outer_workspace_setup_ms,"
               "problem_setup_ms,rhs_setup_ms,outer_min,outer_max,krylov_min,"
               "krylov_max,krylov_sum,dispatched_steps,fine_schur_calls,"
               "preconditioner_calls,coarse_spmv_calls,coarse_iterations,worst_"
               "coarse_true_residual,max_true_residual,solution_error,rhs_per_"
               "second,peak_rss_mib\n";
        trace << "extent,mass,tile_rhs,offset,backend,trial,outer_step,max_"
                 "true_residual,inner_min_iterations,inner_max_iterations,"
                 "dispatched_steps\n";
      }
      int case_index = 0;
      for (auto mass : o.masses)
        for (auto l : o.sizes) {
          int index = case_index++;
          TwoLevel h(l, {o.block, o.block, o.block, o.block}, o.nvec, mass, .25,
                     1, true, o.workers, 8, 1, .65, 1, 1, false);
          SparseCoarse ac(h, mass, .25, 1, true, o.workers);
          auto t = Clock::now();
          CloverSystem c(l, o.tile, mass, .25, 1, true, o.workers, o.self);
          std::vector<Z>().swap(c.w.a.p.x);
          std::vector<Z>().swap(c.parity.eo.p.x);
          std::vector<Z>().swap(c.parity.oe.p.x);
          double problem_setup = seconds(t) * 1000;
          std::cout << "CASE L=" << l << " bare_mass=" << mass
                    << " clover_min_pivot=" << c.min_pivot << std::endl;
          if (o.self || o.validate)
            check_sparse(h, ac, c, metal, o.workers, h.coarse <= 512);
          if (o.validate)
            continue;
          t = Clock::now();
          OuterWorkspace outer(c.w.a, o.workers);
          double outer_setup = seconds(t) * 1000;
          for (int slot = 0; slot < 2; ++slot) {
            bool gpu = (slot + index + int(o.metal_first)) % 2;
            if ((gpu && o.backend == "cpu") || (!gpu && o.backend == "metal"))
              continue;
            for (size_t m = 0; m < methods.size(); ++m) {
              const auto method =
                  methods[(m + index + int(gpu)) % methods.size()];
              @autoreleasepool {
                t = Clock::now();
                if (!gpu) {
                  CPUCloverSchur<float> op(c, o.workers);
                  if (method == "bicgstab") {
                    CPUBiCGStab<float> bi(c.parity.oe, o.workers);
                    uint64_t fine = 0;
                    bi.apply_override = [&](const Cx<float> *x, Cx<float> *y) {
                      ++fine;
                      op.apply(x, y);
                    };
                    HostCloverCorrection<CPUBiCGStab<float>, float> correction(
                        c, bi, o.workers);
                    double setup = seconds(t) * 1000;
                    auto run = [&]() {
                      fine = 0;
                      correction.reset_metrics();
                      auto r = solve(c.w.a, correction, outer, o.workers, 1e-10,
                                     .002, 512);
                      return Outcome{r,
                                     {fine, 0, 0, 0, 0},
                                     correction.prepare_ms,
                                     correction.reconstruct_ms};
                    };
                    compare_measure(c, h, ac, outer, o, method, "cpu", setup,
                                    outer_setup, problem_setup, run, csv,
                                    trace);
                  } else {
                    CPUMultigrid mg(c, h, op, o.workers);
                    mg.coarse_enabled = method == "sparse_mg";
                    std::unique_ptr<CPUCoarseSolve> coarse;
                    if (mg.coarse_enabled) {
                      coarse = std::make_unique<CPUCoarseSolve>(
                          ac, o.tile, o.workers, o.coarse_limit, o.coarse_tol);
                      mg.coarse_override = [&](const Cx<float> *x,
                                               Cx<float> *y) {
                        coarse->solve(x, y);
                      };
                    }
                    CPUFGMRES fg(c, op, mg, o.workers, o.restart);
                    HostCloverCorrection<CPUFGMRES, float> correction(
                        c, fg, o.workers);
                    double setup = seconds(t) * 1000;
                    auto run = [&]() {
                      if (coarse)
                        coarse->reset();
                      mg.calls = mg.schur_calls = fg.fine_calls = 0;
                      correction.reset_metrics();
                      auto r = solve(c.w.a, correction, outer, o.workers, 1e-10,
                                     .002, 256);
                      return Outcome{r,
                                     {mg.schur_calls + fg.fine_calls, mg.calls,
                                      coarse ? coarse->applications : 0,
                                      coarse ? coarse->iterations : 0,
                                      coarse ? coarse->worst_true : 0},
                                     correction.prepare_ms,
                                     correction.reconstruct_ms};
                    };
                    compare_measure(c, h, ac, outer, o, method, "cpu", setup,
                                    outer_setup, problem_setup, run, csv,
                                    trace);
                  }
                } else {
                  GPUCloverSchur op(metal, c);
                  if (method == "bicgstab") {
                    GPUBiCGStab bi(metal, c.parity.oe, o.bi_chunk, o.workers);
                    uint64_t fine = 0;
                    bi.encode_override = [&](id<MTLComputeCommandEncoder> enc,
                                             id<MTLBuffer> x, id<MTLBuffer> y) {
                      ++fine;
                      op.encode(enc, x, y);
                    };
                    DeviceCloverCorrection correction(c, op, bi, o.workers);
                    double setup = seconds(t) * 1000;
                    auto run = [&]() {
                      fine = 0;
                      correction.reset_metrics();
                      auto r = solve(c.w.a, correction, outer, o.workers, 1e-10,
                                     .002, 512);
                      return Outcome{r,
                                     {fine, 0, 0, 0, 0},
                                     correction.prepare_ms,
                                     correction.reconstruct_ms};
                    };
                    compare_measure(c, h, ac, outer, o, method, "metal", setup,
                                    outer_setup, problem_setup, run, csv,
                                    trace);
                  } else {
                    GPUMultigrid mg(c, h, op);
                    mg.coarse_enabled = method == "sparse_mg";
                    std::unique_ptr<GPUCoarseSolve> coarse;
                    if (mg.coarse_enabled) {
                      coarse = std::make_unique<GPUCoarseSolve>(
                          ac, metal, o.tile, o.coarse_limit, o.coarse_tol);
                      mg.coarse_override =
                          [&](id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
                              id<MTLBuffer> y) { coarse->encode(enc, x, y); };
                      mg.coarse_check = [&]() { coarse->check(); };
                    }
                    GPUFGMRES fg(c, op, mg, o.restart);
                    DeviceCloverCorrectionT<GPUFGMRES> correction(c, op, fg,
                                                                  o.workers);
                    double setup = seconds(t) * 1000;
                    auto run = [&]() {
                      if (coarse)
                        coarse->reset();
                      mg.calls = mg.schur_calls = fg.fine_calls = 0;
                      correction.reset_metrics();
                      auto r = solve(c.w.a, correction, outer, o.workers, 1e-10,
                                     .002, 256);
                      return Outcome{r,
                                     {mg.schur_calls + fg.fine_calls, mg.calls,
                                      coarse ? coarse->applications : 0,
                                      coarse ? coarse->iterations : 0,
                                      coarse ? coarse->worst_true : 0},
                                     correction.prepare_ms,
                                     correction.reconstruct_ms};
                    };
                    compare_measure(c, h, ac, outer, o, method, "metal", setup,
                                    outer_setup, problem_setup, run, csv,
                                    trace);
                  }
                }
              }
            }
          }
        }
      if (o.validate)
        std::cout << "PASS: hierarchy and coarse-solver validation; no full "
                     "solves requested."
                  << std::endl;
      else
        std::cout << "PASS: all methods and RHS met full CPU FP64 "
                     "clover-Wilson residual checks."
                  << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
