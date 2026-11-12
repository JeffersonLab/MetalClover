#include "mg_checks.hpp"

using namespace MetalCloverEnv;

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      std::vector<uint32_t> sizes{8}, rhs_values{96, 192};
      std::vector<double> masses{.1, .5};
      uint32_t block = 4, restart = 8, nvec = 6;
      std::array<uint32_t, 4> shape{};
      double gauge = .25, csw = 1, tol = 1e-10, inner_tol = .002, omega = .65,
             setup_omega = .65;
      int workers = 8, samples = 3, sweeps = 8, pre = 1, post = 1;
      bool self = false, apbc = true, validate_only = false, warmup = true;
      std::string csv_path, method_arg;
      for (int j = 1; j < argc; ++j) {
        std::string opt = argv[j];
        auto next = [&]() {
          if (++j >= argc)
            throw std::runtime_error("missing value");
          return std::string(argv[j]);
        };
        if (opt == "--sizes")
          sizes = list_arg(next(), 12);
        else if (opt == "--rhs")
          rhs_values = list_arg(next(), 192);
        else if (opt == "--aggregate")
          block = uint32_t(std::stoul(next()));
        else if (opt == "--aggregate-shape") {
          auto v = list_arg(next(), 12);
          if (v.size() != 4)
            throw std::runtime_error("aggregate-shape needs four dimensions");
          std::copy(v.begin(), v.end(), shape.begin());
        } else if (opt == "--vectors")
          nvec = uint32_t(std::stoul(next()));
        else if (opt == "--setup-omega")
          setup_omega = std::stod(next());
        else if (opt == "--pre-smooth")
          pre = std::stoi(next());
        else if (opt == "--post-smooth")
          post = std::stoi(next());
        else if (opt == "--methods")
          method_arg = next();
        else if (opt == "--no-warmup")
          warmup = false;
        else if (opt == "--mass") {
          masses.clear();
          std::stringstream in(next());
          std::string s;
          while (std::getline(in, s, ','))
            masses.push_back(std::stod(s));
        } else if (opt == "--csw")
          csw = std::stod(next());
        else if (opt == "--gauge-strength")
          gauge = std::stod(next());
        else if (opt == "--tol")
          tol = std::stod(next());
        else if (opt == "--inner-tol")
          inner_tol = std::stod(next());
        else if (opt == "--omega")
          omega = std::stod(next());
        else if (opt == "--setup-sweeps")
          sweeps = std::stoi(next());
        else if (opt == "--restart")
          restart = uint32_t(std::stoul(next()));
        else if (opt == "--samples")
          samples = std::stoi(next());
        else if (opt == "--threads")
          workers = std::stoi(next());
        else if (opt == "--csv")
          csv_path = next();
        else if (opt == "--self-test")
          self = true;
        else if (opt == "--validate-only")
          validate_only = true;
        else if (opt == "--periodic-time")
          apbc = false;
        else if (opt == "--help") {
          std::cout
              << "--sizes 8 --rhs 96,192 --mass .1,.5 --aggregate 4 --restart "
                 "8 --setup-sweeps 8 --omega .65 --samples 3 --threads 8 --csw "
                 "1 --gauge-strength .25 --tol 1e-10 --inner-tol .002 --csv "
                 "results/multigrid.csv --periodic-time --self-test "
                 "--validate-only --vectors 6 --aggregate-shape 4,4,4,4 "
                 "--setup-omega .65 --pre-smooth 1 --post-smooth 1 --methods "
                 "cpu_fgmres_mg,metal_fgmres_mg --no-warmup\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + opt);
      }
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error(
            "set VECLIB_MAXIMUM_THREADS=1 for CPU coarse BLAS");
      if (masses.empty() ||
          !std::all_of(masses.begin(), masses.end(),
                       [](double m) { return m > 0 && std::isfinite(m); }) ||
          !(csw >= 0 && csw <= 3) || !(gauge >= 0 && gauge <= 1) ||
          !(tol >= 1e-13 && tol < 1e-3) ||
          !(inner_tol > tol && inner_tol < .2) ||
          !(omega > 0 && omega <= 1.5) ||
          !(setup_omega > 0 && setup_omega <= 1) || sweeps < 0 || sweeps > 96 ||
          restart < 2 || restart > 12 || workers < 1 || workers > 32 ||
          samples < 1 || samples > 9 || samples % 2 == 0 || block < 2 ||
          block % 2)
        throw std::runtime_error("invalid options");
      if (nvec < 1 || nvec > 24 || pre < 0 || post < 0 || pre + post < 1 ||
          pre + post > 6)
        throw std::runtime_error("invalid basis or smoothing options");
      if (self) {
        sizes = {2, 4};
        rhs_values = {1, 3};
        masses = {.25};
        samples = 1;
        block = 2;
        shape = {};
        restart = 6;
      }
      if (shape == std::array<uint32_t, 4>{})
        shape.fill(block);
      std::string shape_name =
          std::to_string(shape[0]) + "x" + std::to_string(shape[1]) + "x" +
          std::to_string(shape[2]) + "x" + std::to_string(shape[3]);
      for (auto l : sizes) {
        for (auto d : shape)
          if (d < 2 || d % 2 || l % d)
            throw std::runtime_error(
                "even aggregate dimensions must divide even extent");
        for (auto r : rhs_values)
          if (uint64_t(l) * l * l * l * 12 * r > 12000000)
            throw std::runtime_error("FGMRES workspace limit exceeded");
      }
      std::string source;
      for (const char *name : {"kernels.metal", "bicgstab.metal",
                               "clover.metal", "multigrid.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("run from prototype directory");
        source.append((std::istreambuf_iterator<char>(f)), {});
        source += '\n';
      }
      auto t = Clock::now();
      MetalRun metal(source);
      double library_ms = seconds(t) * 1000;
      std::cout << "GPU: " << [[metal.device name] UTF8String] << " OS: " <<
          [[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                << "\nTwo-level odd clover Schur; csw=" << csw
                << " aggregate=" << shape_name << " restart=" << restart
                << " vectors/chirality=" << nvec << " pre/post=" << pre << "/"
                << post << " setup_omega=" << setup_omega
                << " sweeps=" << sweeps << " omega=" << omega
                << " workers=" << workers << " library_ms=" << library_ms
                << std::endl;
      std::ofstream csv, trace;
      if (!csv_path.empty()) {
        csv.open(csv_path);
        trace.open(csv_path + ".trace.csv");
        if (!csv || !trace)
          throw std::runtime_error("cannot create results");
        csv << "extent,rhs,mass,csw,aggregate,coarse_dim,vectors,pre_smooth,"
               "post_smooth,setup_omega,restart,setup_sweeps,omega,backend,"
               "trial,tolerance,total_ms,inner_ms,residual_ms,update_ms,"
               "prepare_ms,reconstruct_ms,hierarchy_setup_ms,workspace_setup_"
               "ms,problem_setup_ms,outer_steps,inner_max_sum,dispatched_steps,"
               "fine_schur_calls,preconditioner_calls,coarse_calls,max_true_"
               "residual,solution_error,ms_per_rhs,rhs_per_second\n";
        trace << "extent,rhs,mass,backend,trial,outer_step,max_true_residual,"
                 "inner_min_iterations,inner_max_iterations,dispatched_steps\n";
      }
      int cases = 0;
      for (auto mass : masses)
        for (auto l : sizes) {
          TwoLevel hierarchy(l, shape, nvec, mass, gauge, csw, apbc, workers,
                             sweeps, omega, setup_omega, pre, post);
          std::cout << "Hierarchy setup L=" << l << " mass=" << mass
                    << " coarse=" << hierarchy.coarse
                    << " ms=" << hierarchy.setup_ms << std::endl;
          for (auto rhs : rhs_values) {
            t = Clock::now();
            CloverSystem c(l, rhs, mass, gauge, csw, apbc, workers, self);
            double problem_ms = seconds(t) * 1000;
            t = Clock::now();
            CPUCloverSchur<float> cpuop(c, workers);
            CPUBiCGStab<float> cb(c.parity.oe, workers);
            uint64_t cpu_bi_calls = 0;
            cb.apply_override = [&](const Cx<float> *x, Cx<float> *y) {
              ++cpu_bi_calls;
              cpuop.apply(x, y);
            };
            HostCloverCorrection<CPUBiCGStab<float>, float> cbd(c, cb, workers);
            double sbc = seconds(t) * 1000;
            t = Clock::now();
            GPUCloverSchur gpuop(metal, c);
            GPUBiCGStab gb(metal, c.parity.oe, 4, workers);
            uint64_t gpu_bi_calls = 0;
            gb.encode_override = [&](id<MTLComputeCommandEncoder> enc,
                                     id<MTLBuffer> x, id<MTLBuffer> y) {
              ++gpu_bi_calls;
              gpuop.encode(enc, x, y);
            };
            DeviceCloverCorrection gbd(c, gpuop, gb, workers);
            double sbg = seconds(t) * 1000;
            t = Clock::now();
            CPUMultigrid cmg(c, hierarchy, cpuop, workers);
            CPUFGMRES cf(c, cpuop, cmg, workers, restart);
            HostCloverCorrection<CPUFGMRES, float> cfd(c, cf, workers);
            double sfc = seconds(t) * 1000;
            t = Clock::now();
            GPUMultigrid gmg(c, hierarchy, gpuop);
            GPUFGMRES gf(c, gpuop, gmg, restart);
            DeviceCloverCorrectionT<GPUFGMRES> gfd(c, gpuop, gf, workers);
            double sfg = seconds(t) * 1000;
            if (self || validate_only)
              check_hierarchy(hierarchy, c, gpuop, gmg, workers);
            if (validate_only) {
              ++cases;
              continue;
            }
            OuterWorkspace outer(c.w.a, workers);
            const std::string names[] = {
                "cpu_bicgstab",        "metal_bicgstab", "cpu_fgmres_smooth",
                "metal_fgmres_smooth", "cpu_fgmres_mg",  "metal_fgmres_mg"};
            double setup[] = {sbc, sbg, sfc, sfg, sfc, sfg};
            std::vector<int> selected;
            if (method_arg.empty())
              selected = {0, 1, 2, 3, 4, 5};
            else {
              std::stringstream in(method_arg);
              std::string name;
              while (std::getline(in, name, ',')) {
                int b = 0;
                while (b < 6 && names[b] != name)
                  ++b;
                if (b == 6)
                  throw std::runtime_error("unknown method " + name);
                if (std::find(selected.begin(), selected.end(), b) !=
                    selected.end())
                  throw std::runtime_error("duplicate method");
                selected.push_back(b);
              }
            }
            auto run = [&](int b) {
              if (b == 0) {
                cpu_bi_calls = 0;
                cbd.reset_metrics();
                return solve(c.w.a, cbd, outer, workers, tol, inner_tol, 512);
              }
              if (b == 1) {
                gpu_bi_calls = 0;
                gbd.reset_metrics();
                return solve(c.w.a, gbd, outer, workers, tol, inner_tol, 512);
              }
              if (b == 2 || b == 4) {
                cmg.coarse_enabled = b == 4;
                cmg.calls = cmg.schur_calls = cf.fine_calls = 0;
                cfd.reset_metrics();
                return solve(c.w.a, cfd, outer, workers, tol, inner_tol, 192);
              }
              gmg.coarse_enabled = b == 5;
              gmg.calls = gmg.schur_calls = gf.fine_calls = 0;
              gfd.reset_metrics();
              return solve(c.w.a, gfd, outer, workers, tol, inner_tol, 192);
            };
            if (!self && warmup)
              for (size_t b = 0; b < selected.size(); ++b)
                run(selected[(b + cases) % selected.size()]);
            for (int trial = 0; trial < samples; ++trial)
              for (size_t b = 0; b < selected.size(); ++b) {
                int method = selected[(b + trial + cases) % selected.size()];
                auto r = run(method);
                double prep = method == 0       ? cbd.prepare_ms
                              : method == 1     ? gbd.prepare_ms
                              : method % 2 == 0 ? cfd.prepare_ms
                                                : gfd.prepare_ms;
                double rec = method == 0       ? cbd.reconstruct_ms
                             : method == 1     ? gbd.reconstruct_ms
                             : method % 2 == 0 ? cfd.reconstruct_ms
                                               : gfd.reconstruct_ms;
                uint64_t fine = method == 0   ? cpu_bi_calls
                                : method == 1 ? gpu_bi_calls
                                : method % 2 == 0
                                    ? cf.fine_calls + cmg.schur_calls
                                    : gf.fine_calls + gmg.schur_calls;
                uint64_t pre_calls = method < 2        ? 0
                                     : method % 2 == 0 ? cmg.calls
                                                       : gmg.calls,
                         coarse = method < 4 ? 0 : pre_calls;
                std::cout << "mass=" << mass << " L=" << l << " RHS=" << rhs
                          << " " << names[method] << " trial=" << trial
                          << " ms=" << r.total_ms << " outer=" << r.outer
                          << " inner=" << r.inner_max_sum << " fineS=" << fine
                          << " Vcycles=" << coarse
                          << " true_r=" << r.max_residual
                          << " x_error=" << r.solution_error << " PASS"
                          << std::endl;
                if (csv) {
                  csv << std::setprecision(12) << l << ',' << rhs << ',' << mass
                      << ',' << csw << ',' << shape_name << ','
                      << hierarchy.coarse << ',' << nvec << ',' << pre << ','
                      << post << ',' << setup_omega << ',' << restart << ','
                      << sweeps << ',' << omega << ',' << names[method] << ','
                      << trial << ',' << tol << ',' << r.total_ms << ','
                      << r.inner_ms << ',' << r.residual_ms << ','
                      << r.update_ms << ',' << prep << ',' << rec << ','
                      << hierarchy.setup_ms << ',' << setup[method] << ','
                      << problem_ms << ',' << r.outer << ',' << r.inner_max_sum
                      << ',' << r.dispatched << ',' << fine << ',' << pre_calls
                      << ',' << coarse << ',' << r.max_residual << ','
                      << r.solution_error << ',' << r.total_ms / rhs << ','
                      << 1000 * rhs / r.total_ms << '\n';
                  csv.flush();
                  for (auto tr : r.trace) {
                    trace << std::setprecision(12) << l << ',' << rhs << ','
                          << mass << ',' << names[method] << ',' << trial << ','
                          << tr.outer << ',' << tr.max_residual << ','
                          << tr.inner.min_iterations << ','
                          << tr.inner.max_iterations << ','
                          << tr.inner.dispatched << '\n';
                  }
                  trace.flush();
                }
              }
            ++cases;
          }
        }
      if (validate_only)
        std::cout << "PASS: " << cases
                  << " hierarchy/operator validation cases." << std::endl;
      else
        std::cout << "PASS: " << cases
                  << " systems; all methods met full CPU FP64 clover-Wilson "
                     "residual checks."
                  << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
