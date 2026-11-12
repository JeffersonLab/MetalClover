#pragma once
#include <cstdint>

namespace MetalCloverEnv {

  // Nominal source-level real operations, FMA=2, complex multiply-add=8.
  // Counts dispatched operator kernels; masked coarse iterations are included.
  // Excludes Krylov vector/reduction kernels, setup SVD, and CPU FP64 checks.
  struct OperatorWork {
    uint64_t fine = 0, coarse = 0, transfer = 0, local_vector = 0;
    uint64_t total() const { return fine + coarse + transfer + local_vector; }
    OperatorWork operator-(const OperatorWork &b) const {
      return {fine - b.fine, coarse - b.coarse, transfer - b.transfer,
              local_vector - b.local_vector};
    }
    double gflops(double elapsed_seconds) const {
      return elapsed_seconds > 0 ? double(total()) / (elapsed_seconds * 1e9)
                                 : 0;
    }
  };

  struct CloverSolveInfo {
    int outer_iterations = 0, krylov_iterations = 0;
    uint64_t fine_schur_applications = 0;
    double seconds = 0, residual_seconds = 0, relative_residual = 0;
    OperatorWork operator_work;
  };

} // namespace MetalCloverEnv
