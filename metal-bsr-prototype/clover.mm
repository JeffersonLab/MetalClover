#include "clover_checks.hpp"
int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      std::vector<uint32_t> sizes{8}, rhs_values{96, 192};
      std::vector<double> masses{.1, .5};
      double csw = 1, gauge = .25, tol = 1e-10, inner_tol = .002;
      int workers = 8, chunk = 4, samples = 3, limit = 512;
      bool self = false, apbc = true;
      std::string csv_path;
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
          rhs_values = list_arg(next(), 384);
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
        else if (opt == "--threads")
          workers = std::stoi(next());
        else if (opt == "--chunk")
          chunk = std::stoi(next());
        else if (opt == "--samples")
          samples = std::stoi(next());
        else if (opt == "--csv")
          csv_path = next();
        else if (opt == "--self-test")
          self = true;
        else if (opt == "--periodic-time")
          apbc = false;
        else if (opt == "--help") {
          std::cout << "--sizes 8 --rhs 96,192 --mass .1,.5 --csw 1 "
                       "--gauge-strength .25 --tol 1e-10 --inner-tol .002 "
                       "--threads 8 --chunk 4 --samples 3 --csv "
                       "results/clover.csv --periodic-time --self-test\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + opt);
      }
      if (masses.empty() ||
          !std::all_of(masses.begin(), masses.end(),
                       [](double m) { return m > 0 && std::isfinite(m); }) ||
          !(csw >= 0 && csw <= 3) || !(gauge >= 0 && gauge <= 1) ||
          !(tol >= 1e-13 && tol < 1e-3) ||
          !(inner_tol > tol && inner_tol < .2) || workers < 1 || workers > 32 ||
          chunk < 1 || chunk > 32 || samples < 1 || samples > 15 ||
          samples % 2 == 0)
        throw std::runtime_error("invalid options");
      if (self) {
        sizes = {2, 4};
        rhs_values = {1, 3, 7};
        masses = {.25};
        samples = 1;
      }
      for (auto l : sizes) {
        if (l % 2)
          throw std::runtime_error("even/odd needs even extents");
        for (auto r : rhs_values)
          if (uint64_t(l) * l * l * l * 12 * r > 18000000)
            throw std::runtime_error("working-set limit exceeded");
      }
      std::string source;
      for (const char *name :
           {"kernels.metal", "bicgstab.metal", "clover.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("run from prototype directory");
        source.append((std::istreambuf_iterator<char>(f)), {});
        source += '\n';
      }
      auto start = Clock::now();
      MetalRun metal(source);
      double library_ms = seconds(start) * 1000;
      std::cout << "GPU: " << [[metal.device name] UTF8String] << " OS: " <<
          [[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                << "\nQDP-basis clover c_sw=" << csw
                << " odd Schur, workers=" << workers << " APBC=" << apbc
                << " tol=" << tol << " inner_tol=" << inner_tol
                << " chunk=" << chunk << " library_ms=" << library_ms
                << std::endl;
      std::ofstream csv, trace;
      if (!csv_path.empty()) {
        csv.open(csv_path);
        trace.open(csv_path + ".trace.csv");
        if (!csv || !trace)
          throw std::runtime_error("cannot create results");
        csv << "extent,rows,rhs,mass,csw,gauge_strength,antiperiodic,backend,"
               "trial,workers,chunk,tolerance,inner_tolerance,total_ms,inner_"
               "ms,residual_ms,update_ms,prepare_ms,reconstruct_ms,workspace_"
               "setup_ms,problem_setup_ms,library_ms,outer_steps,inner_max_sum,"
               "dispatched_steps,max_true_residual,solution_error,ms_per_rhs,"
               "rhs_per_second\n";
        trace << "extent,rhs,mass,backend,trial,outer_step,max_true_residual,"
                 "inner_min_iterations,inner_max_iterations,dispatched_steps\n";
      }
      int cases = 0;
      auto run_case = [&](uint32_t l, uint32_t rhs, double mass,
                          double strength, double coefficient, bool boundary) {
        start = Clock::now();
        CloverSystem c(l, rhs, mass, strength, coefficient, boundary, workers,
                       self);
        double model_ms = seconds(start) * 1000;
        if (self)
          check_clover(c, metal, workers);
        auto &a = c.w.a;
        start = Clock::now();
        CPUCloverSchur<double> op64(c, workers);
        CPUBiCGStab<double> bi64(c.parity.oe, workers);
        bi64.apply_override = [&](const D *x, D *y) { op64.apply(x, y); };
        HostCloverCorrection<CPUBiCGStab<double>, double> cpu64(c, bi64,
                                                                workers);
        double s64 = seconds(start) * 1000;
        start = Clock::now();
        CPUCloverSchur<float> op32(c, workers);
        CPUBiCGStab<float> bi32(c.parity.oe, workers);
        bi32.apply_override = [&](const Cx<float> *x, Cx<float> *y) {
          op32.apply(x, y);
        };
        double s32base = seconds(start) * 1000;
        start = Clock::now();
        HostCloverCorrection<CPUBiCGStab<float>, double> cpu_host(c, bi32,
                                                                  workers);
        double s32h = s32base + seconds(start) * 1000;
        start = Clock::now();
        HostCloverCorrection<CPUBiCGStab<float>, float> cpu_float(c, bi32,
                                                                  workers);
        double s32f = s32base + seconds(start) * 1000;
        start = Clock::now();
        GPUCloverSchur opg(metal, c);
        GPUBiCGStab big(metal, c.parity.oe, chunk, workers);
        big.encode_override = [&](id<MTLComputeCommandEncoder> enc,
                                  id<MTLBuffer> x,
                                  id<MTLBuffer> y) { opg.encode(enc, x, y); };
        double sgbase = seconds(start) * 1000;
        if (self)
          check_device_elimination(c, opg, big, workers);
        start = Clock::now();
        HostCloverCorrection<GPUBiCGStab, double> gpu_host(c, big, workers);
        double sgh = sgbase + seconds(start) * 1000;
        start = Clock::now();
        DeviceCloverCorrection gpu_float(c, opg, big, workers);
        double sgf = sgbase + seconds(start) * 1000;
        start = Clock::now();
        OuterWorkspace outer(a, workers);
        double common_setup = seconds(start) * 1000;
        const std::string names[] = {"cpu_fp64", "cpu_mixed_host64",
                                     "cpu_mixed_fp32", "metal_mixed_host64",
                                     "metal_mixed_fp32"};
        double setup[] = {s64, s32h, s32f, sgh, sgf};
        auto run = [&](int b) {
          if (b == 0) {
            cpu64.reset_metrics();
            return solve(a, cpu64, outer, workers, tol, tol, limit);
          }
          if (b == 1) {
            cpu_host.reset_metrics();
            return solve(a, cpu_host, outer, workers, tol, inner_tol, limit);
          }
          if (b == 2) {
            cpu_float.reset_metrics();
            return solve(a, cpu_float, outer, workers, tol, inner_tol, limit);
          }
          if (b == 3) {
            gpu_host.reset_metrics();
            return solve(a, gpu_host, outer, workers, tol, inner_tol, limit);
          }
          gpu_float.reset_metrics();
          return solve(a, gpu_float, outer, workers, tol, inner_tol, limit);
        };
        if (!self)
          for (int b = 0; b < 5; ++b)
            run((b + cases) % 5);
        for (int trial = 0; trial < samples; ++trial)
          for (int b = 0; b < 5; ++b) {
            int backend = (b + trial + cases) % 5;
            auto r = run(backend);
            double prep[] = {cpu64.prepare_ms, cpu_host.prepare_ms,
                             cpu_float.prepare_ms, gpu_host.prepare_ms,
                             gpu_float.prepare_ms};
            double rec[] = {cpu64.reconstruct_ms, cpu_host.reconstruct_ms,
                            cpu_float.reconstruct_ms, gpu_host.reconstruct_ms,
                            gpu_float.reconstruct_ms};
            std::cout << "mass=" << mass << " L=" << l << " RHS=" << rhs
                      << " csw=" << coefficient << " " << names[backend]
                      << " trial=" << trial << " ms=" << r.total_ms
                      << " outer=" << r.outer << " inner=" << r.inner_max_sum
                      << " dispatched=" << r.dispatched
                      << " true_r=" << r.max_residual
                      << " x_error=" << r.solution_error << " PASS"
                      << std::endl;
            if (csv) {
              csv << std::setprecision(12) << l << ',' << a.p.rows << ',' << rhs
                  << ',' << mass << ',' << coefficient << ',' << strength << ','
                  << boundary << ',' << names[backend] << ',' << trial << ','
                  << workers << ',' << (backend >= 3 ? chunk : 1) << ',' << tol
                  << ',' << (backend == 0 ? tol : inner_tol) << ','
                  << r.total_ms << ',' << r.inner_ms << ',' << r.residual_ms
                  << ',' << r.update_ms << ',' << prep[backend] << ','
                  << rec[backend] << ',' << setup[backend] + common_setup << ','
                  << model_ms << ',' << library_ms << ',' << r.outer << ','
                  << r.inner_max_sum << ',' << r.dispatched << ','
                  << r.max_residual << ',' << r.solution_error << ','
                  << r.total_ms / rhs << ',' << 1000 * rhs / r.total_ms << '\n';
              csv.flush();
              for (auto tr : r.trace) {
                trace << std::setprecision(12) << l << ',' << rhs << ',' << mass
                      << ',' << names[backend] << ',' << trial << ','
                      << tr.outer << ',' << tr.max_residual << ','
                      << tr.inner.min_iterations << ','
                      << tr.inner.max_iterations << ',' << tr.inner.dispatched
                      << '\n';
              }
              trace.flush();
            }
          }
        ++cases;
      };
      for (auto mass : masses)
        for (auto l : sizes)
          for (auto rhs : rhs_values)
            run_case(l, rhs, mass, gauge, csw, apbc);
      if (self) {
        run_case(2, 3, .25, 0, 1, !apbc);
        run_case(2, 3, .25, .25, 0, apbc);
      }
      std::cout << "PASS: " << cases
                << " clover systems; every RHS checked on the full FP64 "
                   "clover-Wilson equation after even reconstruction."
                << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
