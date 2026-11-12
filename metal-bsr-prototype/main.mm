#include "bsr_common.hpp"

using namespace MetalCloverEnv;

int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      std::vector<uint32_t> sizes{4, 8, 12, 16}, rhs_values{1, 4, 12, 48};
      std::string csv_path, mode = "both", shader_path = "kernels.metal";
      int workers = 8, samples = 5;
      double target_ms = 20;
      bool self = false, full = false;
      for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() {
          if (++i >= argc)
            throw std::runtime_error("missing argument");
          return std::string(argv[i]);
        };
        if (a == "--sizes")
          sizes = parse_list(next());
        else if (a == "--rhs")
          rhs_values = parse_list(next());
        else if (a == "--threads")
          workers = std::stoi(next());
        else if (a == "--samples")
          samples = std::stoi(next());
        else if (a == "--target-ms")
          target_ms = std::stod(next());
        else if (a == "--csv")
          csv_path = next();
        else if (a == "--mode")
          mode = next();
        else if (a == "--shader")
          shader_path = next();
        else if (a == "--self-test")
          self = true;
        else if (a == "--full-check")
          full = true;
        else if (a == "--help") {
          std::cout << "Run from prototype directory. Options: --sizes "
                       "4,8,12,16 --rhs 1,4,12,48 --threads 8 --samples 5 "
                       "--target-ms 20 --mode both|dense12|color3_spin4 --csv "
                       "PATH --shader PATH --self-test --full-check\n";
          return 0;
        } else
          throw std::runtime_error("unknown option " + a);
      }
      if (workers < 1 || workers > 32 || samples < 3 || samples % 2 == 0 ||
          !std::isfinite(target_ms) || target_ms <= 0 || target_ms > 1000)
        throw std::runtime_error(
            "require threads 1..32, odd samples >=3, target-ms (0,1000]");
      if (mode != "both" && mode != "dense12" && mode != "color3_spin4")
        throw std::runtime_error("unknown mode");
      // Prevent unexpectedly huge allocations / 32-bit shader index overflow.
      for (auto l : sizes)
        for (auto r : rhs_values)
          if (uint64_t(l) * l * l * l * 12 * r > 100000000)
            throw std::runtime_error("case exceeds prototype memory limit");
      if (!std::getenv("VECLIB_MAXIMUM_THREADS") ||
          std::string(std::getenv("VECLIB_MAXIMUM_THREADS")) != "1")
        throw std::runtime_error(
            "set VECLIB_MAXIMUM_THREADS=1 to avoid nested BLAS threads");
      std::ifstream shader(shader_path);
      if (!shader)
        throw std::runtime_error("cannot open " + shader_path);
      std::string source((std::istreambuf_iterator<char>(shader)), {});
      auto init = Clock::now();
      MetalRun metal(source);
      double init_ms = seconds(init) * 1000;
      std::cout << "GPU: " << [[metal.device name] UTF8String] << " | macOS " <<
          [[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                << "\nCPU workers: " << workers
                << "; Accelerate internal threads: 1; samples: " << samples
                << "; initial Metal library setup: " << init_ms << " ms\n";
      if (self) {
        sizes = {1, 2, 3};
        rhs_values = {1, 3, 4, 7, 12};
        samples = 3;
        full = true;
      }
      std::ofstream csv;
      if (!csv_path.empty()) {
        csv.open(csv_path);
        if (!csv)
          throw std::runtime_error("cannot write CSV");
        csv << "mode,extent,rows,rhs,nnz_blocks,working_set_mib,backend,cpu_"
               "workers,samples,repeats,median_ms,min_ms,max_ms,device_ms,"
               "gflops,checked_outputs,max_abs_error,relative_l2_error,upload_"
               "ms,pipeline_ms\n";
      }
      size_t cases = 0;
      for (bool kron : {false, true}) {
        if (mode != "both" && mode != (kron ? "color3_spin4" : "dense12"))
          continue;
        for (auto extent : sizes)
          for (auto rhs : rhs_values) {
            // Self-tests deliberately include empty rows and missing edges.
            Problem p = make_problem(extent, rhs, kron, self);
            GPUProblem gpu(metal, p);
            std::vector<Z> y(p.outputs());
            std::vector<std::vector<Z>> scratch(
                workers, std::vector<Z>(p.colors * p.spins * p.rhs));
            auto indices = check_indices(p, full);
            std::vector<std::complex<double>> ref;
            ref.reserve(indices.size());
            for (auto idx : indices)
              ref.push_back(reference_one(p, idx));
            struct Result {
              std::string name;
              Timing time;
              Check check;
            };
            std::vector<Result> results;
            auto cpu = [&](std::string name, bool blas, int nth) {
              auto f = [&]() {
                if (blas)
                  cpu_blas(p, y.data(), nth, scratch);
                else
                  cpu_neon(p, y.data(), nth);
              };
              std::fill(y.begin(), y.end(), Z{NAN, NAN});
              f();
              Check check = validate(y.data(), indices, ref);
              Timing t = self ? Timing{} : time_cpu(f, samples, target_ms);
              results.push_back({name, t, check});
            };
            // Alternate CPU/GPU order across cases to reduce systematic order
            // bias.
            auto gpu_result = [&]() {
              std::fill(gpu.output(), gpu.output() + p.outputs(), Z{NAN, NAN});
              gpu.run(1);
              Check check = validate(gpu.output(), indices, ref);
              if (self) {
                results.push_back({"metal", {}, check});
                return;
              }
              results.push_back({"metal_latency",
                                 time_gpu(gpu, samples, target_ms, true),
                                 check});
              results.push_back({"metal_batch",
                                 time_gpu(gpu, samples, target_ms, false),
                                 check});
              // Catch repeated-dispatch or buffer lifetime errors as well.
              validate(gpu.output(), indices, ref);
            };
            if (cases % 2)
              gpu_result();
            cpu("accelerate_1", true, 1);
            cpu("accelerate_parallel", true, workers);
            cpu("neon_parallel", false, workers);
            if (!(cases % 2))
              gpu_result();
            std::cout << p.mode() << " L=" << extent << " rows=" << p.rows
                      << " rhs=" << rhs << " check=" << indices.size() << "/"
                      << p.outputs() << " " << std::flush;
            for (auto &r : results) {
              std::cout << r.name << "=" << std::fixed << std::setprecision(3)
                        << r.time.median_ms << "ms ";
              if (csv)
                csv << std::setprecision(10) << p.mode() << ',' << extent << ','
                    << p.rows << ',' << rhs << ',' << p.col.size() << ','
                    << p.bytes() / 1048576.0 << ',' << r.name << ','
                    << (r.name == "accelerate_1"
                            ? 1
                            : (r.name.find("metal") == 0 ? 0 : workers))
                    << ',' << samples << ',' << r.time.repeats << ','
                    << r.time.median_ms << ',' << r.time.min_ms << ','
                    << r.time.max_ms << ',' << r.time.device_ms << ','
                    << (r.time.median_ms > 0
                            ? p.flops() / (r.time.median_ms * 1e6)
                            : 0)
                    << ',' << r.check.count << ',' << r.check.max_abs << ','
                    << r.check.rel_l2 << ',' << gpu.upload_ms << ','
                    << gpu.compile_ms << '\n';
            }
            std::cout << "PASS\n";
            if (csv)
              csv.flush();
            ++cases;
          }
      }
      std::cout << "PASS: " << cases << " cases; "
                << (full ? "full" : "full small / sampled large")
                << " double-reference validation.\n";
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << '\n';
      return 1;
    }
  }
}
