#include "wilson_common.hpp"

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      std::vector<uint32_t> sizes{8}, rhs_values{96, 192};
      std::vector<double> masses{0.1, 0.5};
      double strength = .25, tol = 1e-10, inner_tol = .002;
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
        } else if (opt == "--gauge-strength")
          strength = std::stod(next());
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
          std::cout
              << "--sizes 8 --rhs 96,192 --mass 0.1,0.5 --gauge-strength .25 "
                 "--tol 1e-10 --inner-tol .002 --threads 8 --chunk 4 --samples "
                 "3 --csv results/wilson.csv --periodic-time --self-test\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + opt);
      }
      if (masses.empty() ||
          !std::all_of(masses.begin(), masses.end(),
                       [](double m) { return m > 0 && std::isfinite(m); }) ||
          !(strength >= 0 && strength <= 1) || !(tol >= 1e-13 && tol < 1e-3) ||
          !(inner_tol > tol && inner_tol < .2) || workers < 1 || workers > 32 ||
          chunk < 1 || chunk > 32 || samples < 1 || samples > 15 ||
          samples % 2 == 0)
        throw std::runtime_error("invalid options");
      if (self) {
        sizes = {1, 2, 3};
        rhs_values = {1, 3, 7};
        masses = {.25};
        samples = 1;
      }
      for (auto l : sizes)
        for (auto r : rhs_values)
          if (uint64_t(l) * l * l * l * 12 * r > 18000000)
            throw std::runtime_error("working-set limit exceeded");
      std::ifstream k("kernels.metal"), s("bicgstab.metal");
      if (!k || !s)
        throw std::runtime_error("run from prototype directory");
      std::string source((std::istreambuf_iterator<char>(k)), {});
      source += '\n';
      source.append((std::istreambuf_iterator<char>(s)), {});
      auto init = Clock::now();
      MetalRun metal(source);
      double library_ms = seconds(init) * 1000;
      std::cout << "GPU: " << [[metal.device name] UTF8String] << " OS: " <<
          [[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                << "\nWilson BiCGStab workers=" << workers
                << " library_ms=" << library_ms << " APBC=" << apbc
                << " gauge=" << strength << " tol=" << tol
                << " inner_tol=" << inner_tol << " chunk=" << chunk
                << std::endl;
      std::ofstream csv, trace;
      if (!csv_path.empty()) {
        csv.open(csv_path);
        trace.open(csv_path + ".trace.csv");
        if (!csv || !trace)
          throw std::runtime_error("cannot create results");
        csv << "extent,rows,rhs,mass,gauge_strength,antiperiodic,backend,trial,"
               "workers,chunk,tolerance,inner_tolerance,total_ms,inner_ms,"
               "residual_ms,update_ms,setup_ms,library_ms,outer_steps,inner_"
               "max_sum,dispatched_steps,max_true_residual,solution_error,ms_"
               "per_rhs,rhs_per_second\n";
        trace << "extent,rhs,mass,backend,trial,outer_step,max_true_residual,"
                 "inner_min_iterations,inner_max_iterations,dispatched_steps\n";
      }
      int cases = 0;
      for (auto mass : masses)
        for (auto l : sizes)
          for (auto rhs : rhs_values) {
            auto w = make_wilson(l, rhs, mass, strength, apbc, workers, self);
            auto &a = w.a;
            if (self) {
              validate_wilson(w, metal, workers);
              auto free = make_wilson(l, rhs, mass, 0, apbc, workers);
              validate_wilson(free, metal, workers);
            }
            auto t = Clock::now();
            CPUBiCGStab<float> cf(a, workers);
            double sf = seconds(t) * 1000;
            t = Clock::now();
            CPUBiCGStab<double> cd(a, workers);
            double sd = seconds(t) * 1000;
            t = Clock::now();
            GPUBiCGStab gpu(metal, a, chunk, workers);
            double sg = seconds(t) * 1000;
            t = Clock::now();
            OuterWorkspace outer(a, workers);
            double common_setup = seconds(t) * 1000;
            auto run = [&](int b) {
              if (b == 0)
                return solve(a, cd, outer, workers, tol, tol, limit);
              if (b == 1)
                return solve(a, cf, outer, workers, tol, inner_tol, limit);
              return solve(a, gpu, outer, workers, tol, inner_tol, limit);
            };
            std::string names[] = {"cpu_fp64", "cpu_mixed", "metal_mixed"};
            double setup[] = {sd + common_setup, sf + common_setup,
                              sg + common_setup};
            if (!self)
              for (int b = 0; b < 3; ++b)
                run((b + cases) % 3);
            for (int trial = 0; trial < samples; ++trial)
              for (int b = 0; b < 3; ++b) {
                int backend = (b + trial + cases) % 3;
                auto r = run(backend);
                std::cout << "mass=" << mass << " L=" << l << " RHS=" << rhs
                          << " " << names[backend] << " trial=" << trial
                          << " ms=" << r.total_ms << " outer=" << r.outer
                          << " inner=" << r.inner_max_sum
                          << " dispatched=" << r.dispatched
                          << " true_r=" << r.max_residual
                          << " x_error=" << r.solution_error << " PASS"
                          << std::endl;
                if (csv) {
                  csv << std::setprecision(12) << l << ',' << a.p.rows << ','
                      << rhs << ',' << mass << ',' << strength << ',' << apbc
                      << ',' << names[backend] << ',' << trial << ',' << workers
                      << ',' << (backend == 2 ? chunk : 1) << ',' << tol << ','
                      << (backend == 0 ? tol : inner_tol) << ',' << r.total_ms
                      << ',' << r.inner_ms << ',' << r.residual_ms << ','
                      << r.update_ms << ',' << setup[backend] << ','
                      << library_ms << ',' << r.outer << ',' << r.inner_max_sum
                      << ',' << r.dispatched << ',' << r.max_residual << ','
                      << r.solution_error << ',' << r.total_ms / rhs << ','
                      << 1000 * rhs / r.total_ms << '\n';
                  csv.flush();
                  for (auto tr : r.trace) {
                    trace << std::setprecision(12) << l << ',' << rhs << ','
                          << mass << ',' << names[backend] << ',' << trial
                          << ',' << tr.outer << ',' << tr.max_residual << ','
                          << tr.inner.min_iterations << ','
                          << tr.inner.max_iterations << ','
                          << tr.inner.dispatched << '\n';
                  }
                  trace.flush();
                }
              }
            ++cases;
          }
      std::cout << "PASS: " << cases
                << " Wilson systems; full true FP64 residual and "
                   "known-solution checks on every RHS."
                << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
