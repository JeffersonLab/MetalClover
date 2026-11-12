#pragma once
#include "fgmres.hpp"

namespace MetalCloverEnv {

  void check_hierarchy(const TwoLevel &h, const CloverSystem &c,
                       GPUCloverSchur &op, GPUMultigrid &gm, int workers) {
    double err = 0;
    for (uint32_t a = 0; a < h.aggregates; ++a)
      for (uint32_t ch = 0; ch < 2; ++ch)
        for (uint32_t i = 0; i < h.nvec; ++i)
          for (uint32_t j = 0; j < h.nvec; ++j) {
            CD sum = 0;
            for (uint32_t p = h.offsets[a]; p < h.offsets[a + 1]; ++p)
              for (uint32_t color = 0; color < 3; ++color)
                for (uint32_t s = 2 * ch; s < 2 * ch + 2; ++s) {
                  size_t q = (size_t(h.sites[p]) * 12 + color * 4 + s) * h.nvec;
                  auto x = h.basis[q + i], y = h.basis[q + j];
                  sum += std::conj(CD(x.re, x.im)) * CD(y.re, y.im);
                }
            err = std::max(err, std::abs(sum - CD(i == j ? 1 : 0)));
          }
    if (err > 2e-12)
      throw std::runtime_error("P^H P != I");
    for (uint32_t i = 0; i < h.coarse; ++i)
      for (uint32_t j = 0; j < h.coarse; ++j) {
        auto a = h.ac64[size_t(i) * h.coarse + j],
             b = h.ac64[size_t(j) * h.coarse + i];
        double sign =
            (i % (2 * h.nvec) / h.nvec) == (j % (2 * h.nvec) / h.nvec) ? 1 : -1;
        if (std::abs(CD(a.re, a.im) - sign * std::conj(CD(b.re, b.im))) > 2e-11)
          throw std::runtime_error("coarse gamma5 Hermiticity failed");
      }
    uint32_t rhs = c.w.a.p.rhs;
    std::vector<D> yc(size_t(h.coarse) * rhs), fine(c.parity.oe.p.outputs()),
        sf(fine.size()), coarse(yc.size()), exact(yc.size()),
        recovered(yc.size());
    std::mt19937_64 rng(88291);
    std::uniform_real_distribution<double> random(-1, 1);
    for (auto &z : yc)
      z = {random(rng), random(rng)};
    h.prolong(yc.data(), fine.data(), rhs);
    CPUCloverSchur<double> fineop(c, workers);
    fineop.apply(fine.data(), sf.data());
    h.restrict_vector(sf.data(), coarse.data(), rhs);
    for (uint32_t i = 0; i < h.coarse; ++i)
      for (uint32_t r = 0; r < rhs; ++r) {
        CD sum = 0;
        for (uint32_t j = 0; j < h.coarse; ++j) {
          auto a = h.ac64[size_t(i) * h.coarse + j],
               b = yc[size_t(j) * rhs + r];
          sum += CD(a.re, a.im) * CD(b.re, b.im);
        }
        exact[size_t(i) * rhs + r] = {sum.real(), sum.imag()};
      }
    err = 0;
    for (size_t j = 0; j < coarse.size(); ++j)
      err = std::max(err, std::hypot(coarse[j].re - exact[j].re,
                                     coarse[j].im - exact[j].im));
    if (err > 3e-11)
      throw std::runtime_error("Galerkin Ac y != P^H S P y");
    for (uint32_t i = 0; i < h.coarse; ++i)
      for (uint32_t r = 0; r < rhs; ++r) {
        CD sum = 0;
        for (uint32_t j = 0; j < h.coarse; ++j) {
          auto a = h.inverse64[size_t(i) * h.coarse + j],
               b = coarse[size_t(j) * rhs + r];
          sum += CD(a.re, a.im) * CD(b.re, b.im);
        }
        recovered[size_t(i) * rhs + r] = {sum.real(), sum.imag()};
      }
    err = 0;
    for (size_t j = 0; j < yc.size(); ++j)
      err = std::max(err, std::hypot(yc[j].re - recovered[j].re,
                                     yc[j].im - recovered[j].im));
    if (err > 5e-10)
      throw std::runtime_error("coarse correction failed on range(SP)");
    auto xf = block_cast<float>(sf);
    std::vector<Cx<float>> cr(size_t(h.coarse) * rhs), py(fine.size());
    h.restrict_vector(xf.data(), cr.data(), rhs);
    h.prolong(cr.data(), py.data(), rhs);
    auto in = op.oe.buffers[3], out = op.oe.buffers[4];
    std::memcpy([in contents], xf.data(), xf.size() * 8);
    @autoreleasepool {
      auto cmd = [op.m.queue commandBuffer];
      auto enc = [cmd computeCommandEncoder];
      gm.dispatch(enc, 1, {in, gm.basis, gm.offsets, gm.sites, gm.coarse_rhs},
                  size_t(h.coarse) * rhs);
      gm.dispatch(enc, 3, {gm.coarse_rhs, gm.basis, gm.aggregate, out},
                  fine.size());
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
      if (cmd.status == MTLCommandBufferStatusError)
        throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);
    }
    auto *gp = static_cast<Cx<float> *>([out contents]);
    auto *gr = static_cast<Cx<float> *>([gm.coarse_rhs contents]);
    err = 0;
    for (size_t j = 0; j < cr.size(); ++j)
      err = std::max(
          err, double(std::hypot(cr[j].re - gr[j].re, cr[j].im - gr[j].im)));
    for (size_t j = 0; j < py.size(); ++j)
      err = std::max(
          err, double(std::hypot(py[j].re - gp[j].re, py[j].im - gp[j].im)));
    if (err > 2e-4)
      throw std::runtime_error(
          "CPU/Metal restriction or prolongation mismatch");
    CPUCloverSchur<float> cpuop(c, workers);
    CPUMultigrid cmg(c, h, cpuop, workers);
    std::vector<Cx<float>> cz(fine.size());
    cmg.apply(xf.data(), cz.data());
    @autoreleasepool {
      auto cmd = [op.m.queue commandBuffer];
      auto enc = [cmd computeCommandEncoder];
      gm.encode(enc, in, out);
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
      if (cmd.status == MTLCommandBufferStatusError)
        throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);
    }
    double e2 = 0, n2 = 0;
    for (size_t j = 0; j < cz.size(); ++j) {
      e2 += std::pow(cz[j].re - gp[j].re, 2) + std::pow(cz[j].im - gp[j].im, 2);
      n2 += cz[j].re * cz[j].re + cz[j].im * cz[j].im;
    }
    if (std::sqrt(e2 / std::max(n2, 1e-300)) > 2e-5)
      throw std::runtime_error("CPU/Metal V-cycle mismatch");
    std::cout << "MG orthogonality/gamma5/Galerkin/coarse correction/CPU-Metal "
                 "transfers/V-cycle checks PASS; coarse="
              << h.coarse << std::endl;
  }

} // namespace MetalCloverEnv
