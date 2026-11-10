#pragma once
#include "sparse_coarse.hpp"

void check_sparse(const TwoLevel &h, const SparseCoarse &a,
                  const CloverSystem &c, MetalRun &metal, int workers,
                  bool dense_oracle) {
  uint32_t rhs = c.w.a.p.rhs;
  std::mt19937_64 rng(71426);
  std::uniform_real_distribution<double> rnd(-1, 1);
  std::vector<D> yc(size_t(a.n) * rhs), fine(c.parity.oe.p.outputs()),
      sf(fine.size()), rp(yc.size()), ac(yc.size());
  for (auto &z : yc)
    z = {rnd(rng), rnd(rng)};
  h.prolong(yc.data(), fine.data(), rhs);
  CPUCloverSchur<double> sop(c, workers);
  sop.apply(fine.data(), sf.data());
  h.restrict_vector(sf.data(), rp.data(), rhs);
  a.apply(yc.data(), ac.data(), rhs);
  double err = 0, norm = 0;
  for (size_t j = 0; j < ac.size(); ++j) {
    err += std::norm(CD(rp[j].re - ac[j].re, rp[j].im - ac[j].im));
    norm += std::norm(CD(rp[j].re, rp[j].im));
  }
  const double galerkin_error = std::sqrt(err / norm);
  if (!std::isfinite(galerkin_error) || galerkin_error > 2e-12)
    throw std::runtime_error("sparse Galerkin != R S P");
  double herm = 0;
  for (uint32_t ar = 0; ar < a.rows; ++ar)
    for (uint32_t e = a.row[ar]; e < a.row[ar + 1]; ++e) {
      uint32_t b = a.col[e];
      auto it = std::lower_bound(a.col.begin() + a.row[b],
                                 a.col.begin() + a.row[b + 1], ar);
      if (it == a.col.begin() + a.row[b + 1] || *it != ar)
        throw std::runtime_error("asymmetric coarse pattern");
      size_t rev = it - a.col.begin();
      for (uint32_t i = 0; i < a.bs; ++i)
        for (uint32_t j = 0; j < a.bs; ++j) {
          auto x = a.a64[(size_t(e) * a.bs + i) * a.bs + j],
               y = a.a64[(rev * a.bs + j) * a.bs + i];
          double sign = (i < h.nvec) == (j < h.nvec) ? 1 : -1;
          herm = std::max(herm, std::abs(CD(x.re, x.im) -
                                         sign * std::conj(CD(y.re, y.im))));
        }
    }
  if (herm > 3e-11)
    throw std::runtime_error("sparse gamma5 Hermiticity failed");
  if (dense_oracle) {
    TwoLevel old(h.extent, h.shape, h.nvec, c.w.a.mass, c.w.strength,
                 c.coefficient, c.w.antiperiodic, workers, h.sweeps, h.omega,
                 h.setup_omega, h.pre_smooth, h.post_smooth);
    std::vector<D> dense(size_t(a.n) * a.n);
    for (uint32_t ar = 0; ar < a.rows; ++ar)
      for (uint32_t e = a.row[ar]; e < a.row[ar + 1]; ++e)
        for (uint32_t i = 0; i < a.bs; ++i)
          for (uint32_t j = 0; j < a.bs; ++j)
            dense[size_t(ar * a.bs + i) * a.n + a.col[e] * a.bs + j] =
                a.a64[(size_t(e) * a.bs + i) * a.bs + j];
    double delta = 0;
    for (size_t j = 0; j < dense.size(); ++j)
      delta = std::max(delta, std::hypot(dense[j].re - old.ac64[j].re,
                                         dense[j].im - old.ac64[j].im));
    if (delta > 2e-11)
      throw std::runtime_error(
          "sparse assembly != independent dense column probing");
  }
  auto xf = block_cast<float>(yc);
  std::vector<Cx<float>> cpu(xf.size());
  a.apply(xf.data(), cpu.data(), rhs);
  GPUCoarseSolve gpu(a, metal, rhs, 64, 1e-5);
  auto in = [metal.device newBufferWithBytes:xf.data()
                                      length:xf.size() * 8
                                     options:MTLResourceStorageModeShared];
  auto out = [metal.device newBufferWithLength:xf.size() * 8
                                       options:MTLResourceStorageModeShared];
  auto finish = [&](id<MTLCommandBuffer> cmd,
                    id<MTLComputeCommandEncoder> enc) {
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
    if (cmd.status == MTLCommandBufferStatusError)
      throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);
  };
  @autoreleasepool {
    auto cmd = [metal.queue commandBuffer];
    auto enc = [cmd computeCommandEncoder];
    gpu.apply(enc, in, out, false);
    finish(cmd, enc);
  }
  auto *y = static_cast<Cx<float> *>([out contents]);
  err = norm = 0;
  for (size_t j = 0; j < cpu.size(); ++j) {
    err += std::norm(CD(cpu[j].re - y[j].re, cpu[j].im - y[j].im));
    norm += std::norm(CD(cpu[j].re, cpu[j].im));
  }
  if (std::sqrt(err / norm) > 3e-6)
    throw std::runtime_error("CPU/Metal sparse multiply mismatch");
  CPUCoarseSolve cb(a, rhs, workers, 64, 1e-5);
  cb.solve(xf.data(), cpu.data());
  @autoreleasepool {
    auto cmd = [metal.queue commandBuffer];
    auto enc = [cmd computeCommandEncoder];
    gpu.encode(enc, in, out);
    finish(cmd, enc);
    gpu.check();
  }
  err = norm = 0;
  for (size_t j = 0; j < cpu.size(); ++j) {
    err += std::norm(CD(cpu[j].re - y[j].re, cpu[j].im - y[j].im));
    norm += std::norm(CD(cpu[j].re, cpu[j].im));
  }
  if (std::sqrt(err / norm) > 1e-4 || cb.worst_true > 1e-3 ||
      gpu.worst_true > 1e-3)
    throw std::runtime_error("coarse iterative CPU/Metal convergence mismatch");
  std::cout << "Sparse Galerkin/gamma5/CPU-Metal SpMV/coarse solve PASS; "
               "dense_oracle="
            << (dense_oracle ? "checked" : "skipped")
            << " relative_Galerkin=" << galerkin_error
            << " coarse_solution_CPU_Metal=" << std::sqrt(err / norm)
            << " coarse_true_cpu=" << cb.worst_true
            << " metal=" << gpu.worst_true << std::endl;
}
