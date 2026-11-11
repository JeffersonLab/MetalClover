#include "clover_checks.hpp"
int main() {
  @autoreleasepool {
    try {
      std::string source;
      for (const char *name :
           {"kernels.metal", "bicgstab.metal", "clover.metal"}) {
        std::ifstream f(name);
        if (!f)
          throw std::runtime_error("missing shader");
        source.append((std::istreambuf_iterator<char>(f)), {});
        source += '\n';
      }
      MetalRun metal(source);
      int workers = 8;
      for (uint32_t l : {2, 4})
        for (uint32_t rhs : {1, 3, 7}) {
          CloverSystem c(l, rhs, .25, .25, 1, true, workers, true);
          check_clover(c, metal, workers);
          CPUCloverSchur<double> op64(c, workers);
          CPUBiCGStab<double> bi64(c.parity.oe, workers);
          bi64.apply_override = [&](const D *x, D *y) { op64.apply(x, y); };
          HostCloverCorrection<CPUBiCGStab<double>, double> cpu(c, bi64,
                                                                workers);
          GPUCloverSchur opg(metal, c);
          GPUBiCGStab big(metal, c.parity.oe, 4, workers);
          big.encode_override = [&](id<MTLComputeCommandEncoder> enc,
                                    id<MTLBuffer> x,
                                    id<MTLBuffer> y) { opg.encode(enc, x, y); };
          HostCloverCorrection<GPUBiCGStab, double> gpu(c, big, workers);
          OuterWorkspace outer(c.w.a, workers);
          auto a = solve(c.w.a, cpu, outer, workers, 1e-10, 1e-10, 512);
          auto b = solve(c.w.a, gpu, outer, workers, 1e-10, .002, 512);
          std::cout << "Clover stage-1 full CPU and Metal/FP64-elimination "
                       "solves PASS residuals="
                    << a.max_residual << "," << b.max_residual << std::endl;
        }
      CloverSystem free(2, 3, .25, 0, 1, true, workers),
          zero(2, 3, .25, .25, 0, true, workers);
      check_clover(free, metal, workers);
      check_clover(zero, metal, workers);
      for (auto *c : {&free, &zero})
        for (auto z : c->w.a.local64)
          if (std::hypot(z.re, z.im) > 1e-14)
            throw std::runtime_error("free/zero-csw clover failed");
      std::cout
          << "PASS: clover stage 1 validated before GPU elimination changes."
          << std::endl;
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
