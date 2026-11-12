#pragma once
#include "clover_checks.hpp"

namespace MetalCloverEnv {

  std::vector<CD> dense_inverse(std::vector<CD> a, uint32_t n) {
    std::vector<CD> b(size_t(n) * n);
    for (uint32_t i = 0; i < n; ++i)
      b[size_t(i) * n + i] = 1;
    for (uint32_t k = 0; k < n; ++k) {
      uint32_t p = k;
      for (uint32_t i = k + 1; i < n; ++i)
        if (std::abs(a[size_t(i) * n + k]) > std::abs(a[size_t(p) * n + k]))
          p = i;
      if (std::abs(a[size_t(p) * n + k]) < 1e-13)
        throw std::runtime_error("singular Galerkin coarse operator");
      for (uint32_t j = 0; j < n; ++j) {
        std::swap(a[size_t(k) * n + j], a[size_t(p) * n + j]);
        std::swap(b[size_t(k) * n + j], b[size_t(p) * n + j]);
      }
      CD pivot = a[size_t(k) * n + k];
      for (uint32_t j = 0; j < n; ++j) {
        a[size_t(k) * n + j] /= pivot;
        b[size_t(k) * n + j] /= pivot;
      }
      for (uint32_t i = 0; i < n; ++i)
        if (i != k) {
          CD z = a[size_t(i) * n + k];
          for (uint32_t j = 0; j < n; ++j) {
            a[size_t(i) * n + j] -= z * a[size_t(k) * n + j];
            b[size_t(i) * n + j] -= z * b[size_t(k) * n + j];
          }
        }
    }
    return b;
  }

  struct TwoLevel {
    uint32_t extent, rows, aggregates, coarse, nvec;
    std::array<uint32_t, 4> shape;
    int workers, sweeps, pre_smooth, post_smooth;
    double omega, setup_omega, setup_ms = 0;
    std::vector<uint32_t> aggregate, offsets, sites;
    std::vector<D> basis, ac64, inverse64;
    std::vector<Cx<float>> basis32, inverse32;
    TwoLevel(uint32_t l, std::array<uint32_t, 4> b, uint32_t nv, double mass,
             double gauge, double csw, bool apbc, int workers, int sweeps,
             double omega, double setup_omega, int pre, int post,
             bool build_dense = true)
        : extent(l), nvec(nv), shape(b), workers(workers), sweeps(sweeps),
          pre_smooth(pre), post_smooth(post), omega(omega),
          setup_omega(setup_omega) {
      auto start = Clock::now();
      aggregates = 1;
      for (auto d : shape) {
        if (d < 2 || d % 2 || l % d)
          throw std::runtime_error(
              "aggregate dimensions must be even and divide extent");
        aggregates *= l / d;
      }
      coarse = aggregates * 2 * nvec;
      if (build_dense && coarse > 512)
        throw std::runtime_error("dense coarse prototype limited to 512 "
                                 "unknowns; increase aggregate "
                                 "size or reduce vectors");
      CloverSystem setup(l, nvec, mass, gauge, csw, apbc, workers);
      rows = setup.parity.oe.p.rows;
      aggregate.resize(rows);
      std::vector<std::vector<uint32_t>> lists(aggregates);
      for (uint32_t i = 0; i < rows; ++i) {
        uint32_t x = setup.parity.odd[i], a = 0, stride = 1;
        for (int mu = 0; mu < 4; ++mu) {
          a += ((x % l) / shape[mu]) * stride;
          x /= l;
          stride *= l / shape[mu];
        }
        aggregate[i] = a;
        lists[a].push_back(i);
      }
      offsets.push_back(0);
      for (auto &list : lists) {
        sites.insert(sites.end(), list.begin(), list.end());
        offsets.push_back(uint32_t(sites.size()));
      }
      std::mt19937_64 rng(771337);
      std::normal_distribution<double> random(0, 1);
      std::vector<CD> constants(12 * nvec);
      for (auto &z : constants)
        z = {random(rng), random(rng)};
      std::vector<D> test(size_t(rows) * 12 * nvec), ap(test.size()),
          relax(test.size()), odd_inverse(setup.oddA.size());
      for (size_t i = 0; i < setup.oddA.size() / 36; ++i) {
        M6 a;
        for (int j = 0; j < 36; ++j) {
          auto z = setup.oddA[i * 36 + j];
          a[j] = {z.re, z.im};
        }
        double pivot = 1e100;
        auto inv = invert_cholesky(a, pivot);
        for (int j = 0; j < 36; ++j)
          odd_inverse[i * 36 + j] = {inv[j].real(), inv[j].imag()};
      }
      for (uint32_t i = 0; i < rows; ++i)
        for (uint32_t q = 0; q < 12; ++q)
          for (uint32_t v = 0; v < nvec; ++v) {
            auto z =
                constants[q * nvec + v] + .2 * CD(random(rng), random(rng));
            test[(size_t(i) * 12 + q) * nvec + v] = {z.real(), z.imag()};
          }
      CPUCloverSchur<double> op(setup, workers);
      for (int it = 0; it < sweeps; ++it) {
        op.apply(test.data(), ap.data());
        local_blocks(odd_inverse, ap.data(), relax.data(), rows, nvec, workers);
        for (size_t j = 0; j < test.size(); ++j) {
          test[j].re -= setup_omega * relax[j].re;
          test[j].im -= setup_omega * relax[j].im;
        }
      }
      basis = test;
      // Twice-reorthogonalized local QR, separately in each gamma5 eigenspace.
      for (uint32_t a = 0; a < aggregates; ++a)
        for (uint32_t ch = 0; ch < 2; ++ch)
          for (uint32_t v = 0; v < nvec; ++v) {
            for (int pass = 0; pass < 2; ++pass)
              for (uint32_t k = 0; k < v; ++k) {
                CD h = 0;
                for (uint32_t p = offsets[a]; p < offsets[a + 1]; ++p)
                  for (uint32_t c = 0; c < 3; ++c)
                    for (uint32_t s = ch * 2; s < ch * 2 + 2; ++s) {
                      size_t q = (size_t(sites[p]) * 12 + c * 4 + s) * nvec;
                      auto x = basis[q + k], y = basis[q + v];
                      h += std::conj(CD(x.re, x.im)) * CD(y.re, y.im);
                    }
                for (uint32_t p = offsets[a]; p < offsets[a + 1]; ++p)
                  for (uint32_t c = 0; c < 3; ++c)
                    for (uint32_t s = ch * 2; s < ch * 2 + 2; ++s) {
                      size_t q = (size_t(sites[p]) * 12 + c * 4 + s) * nvec;
                      auto x = basis[q + k];
                      auto z = h * CD(x.re, x.im);
                      basis[q + v].re -= z.real();
                      basis[q + v].im -= z.imag();
                    }
              }
            double norm = 0;
            for (uint32_t p = offsets[a]; p < offsets[a + 1]; ++p)
              for (uint32_t c = 0; c < 3; ++c)
                for (uint32_t s = ch * 2; s < ch * 2 + 2; ++s) {
                  auto z =
                      basis[(size_t(sites[p]) * 12 + c * 4 + s) * nvec + v];
                  norm += z.re * z.re + z.im * z.im;
                }
            if (norm < 1e-24)
              throw std::runtime_error("rank-deficient local multigrid basis");
            norm = std::sqrt(norm);
            for (uint32_t p = offsets[a]; p < offsets[a + 1]; ++p)
              for (uint32_t c = 0; c < 3; ++c)
                for (uint32_t s = ch * 2; s < ch * 2 + 2; ++s) {
                  auto &z =
                      basis[(size_t(sites[p]) * 12 + c * 4 + s) * nvec + v];
                  z.re /= norm;
                  z.im /= norm;
                }
          }
      basis32 = block_cast<float>(basis);
      if (!build_dense) {
        setup_ms = seconds(start) * 1000;
        return;
      }
      CloverSystem columns(l, coarse, mass, gauge, csw, apbc, workers);
      CPUCloverSchur<double> sop(columns, workers);
      std::vector<D> pmat(size_t(rows) * 12 * coarse), sp(pmat.size());
      for (uint32_t i = 0; i < rows; ++i)
        for (uint32_t q = 0; q < 12; ++q)
          for (uint32_t v = 0; v < nvec; ++v)
            pmat[(size_t(i) * 12 + q) * coarse + aggregate[i] * 2 * nvec +
                 (q % 4 / 2) * nvec + v] =
                basis[(size_t(i) * 12 + q) * nvec + v];
      sop.apply(pmat.data(), sp.data());
      ac64.resize(size_t(coarse) * coarse);
      restrict_vector(sp.data(), ac64.data(), coarse);
      std::vector<CD> ac(ac64.size());
      for (size_t j = 0; j < ac.size(); ++j)
        ac[j] = {ac64[j].re, ac64[j].im};
      auto inv = dense_inverse(ac, coarse);
      inverse64.resize(inv.size());
      for (size_t j = 0; j < inv.size(); ++j)
        inverse64[j] = {inv[j].real(), inv[j].imag()};
      inverse32 = block_cast<float>(inverse64);
      setup_ms = seconds(start) * 1000;
    }
    template <class T> const std::vector<Cx<T>> &weights() const {
      if constexpr (std::is_same_v<T, double>)
        return basis;
      else
        return basis32;
    }
    template <class T>
    void restrict_vector(const Cx<T> *x, Cx<T> *y, uint32_t rhs) const {
      auto &p = weights<T>();
      using V = SIMD<T>;
      chunks(coarse, workers, [&](int, uint32_t first, uint32_t last) {
        for (uint32_t row = first; row < last; ++row) {
          uint32_t a = row / (2 * nvec), ch = row % (2 * nvec) / nvec,
                   v = row % nvec, r = 0;
          for (; r + V::width <= rhs; r += V::width) {
            auto re = V::zero(), im = V::zero();
            for (uint32_t j = offsets[a]; j < offsets[a + 1]; ++j)
              for (uint32_t c = 0; c < 3; ++c)
                for (uint32_t s = 2 * ch; s < 2 * ch + 2; ++s) {
                  size_t q = size_t(sites[j]) * 12 + c * 4 + s;
                  auto w = p[q * nvec + v];
                  auto z = V::load(x + q * rhs + r);
                  re = V::ma(re, z.first, w.re);
                  re = V::ma(re, z.second, w.im);
                  im = V::ma(im, z.second, w.re);
                  im = V::ms(im, z.first, w.im);
                }
            V::store(y + size_t(row) * rhs + r, re, im);
          }
          for (; r < rhs; ++r) {
            Cx<T> sum{};
            for (uint32_t j = offsets[a]; j < offsets[a + 1]; ++j)
              for (uint32_t c = 0; c < 3; ++c)
                for (uint32_t s = 2 * ch; s < 2 * ch + 2; ++s) {
                  size_t q = size_t(sites[j]) * 12 + c * 4 + s;
                  auto w = p[q * nvec + v], z = x[q * rhs + r];
                  sum.re += w.re * z.re + w.im * z.im;
                  sum.im += w.re * z.im - w.im * z.re;
                }
            y[size_t(row) * rhs + r] = sum;
          }
        }
      });
    }
    template <class T>
    void prolong(const Cx<T> *x, Cx<T> *y, uint32_t rhs,
                 bool add = false) const {
      auto &p = weights<T>();
      using V = SIMD<T>;
      chunks(rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
        for (uint32_t q = first; q < last; ++q) {
          uint32_t base = aggregate[q / 12] * 2 * nvec + (q % 4 / 2) * nvec,
                   r = 0;
          for (; r + V::width <= rhs; r += V::width) {
            auto re = V::zero(), im = V::zero();
            for (uint32_t v = 0; v < nvec; ++v) {
              auto w = p[size_t(q) * nvec + v];
              auto z = V::load(x + size_t(base + v) * rhs + r);
              re = V::ma(re, z.first, w.re);
              re = V::ms(re, z.second, w.im);
              im = V::ma(im, z.second, w.re);
              im = V::ma(im, z.first, w.im);
            }
            if (add) {
              auto z = V::load(y + size_t(q) * rhs + r);
              re = V::ma(re, z.first, T(1));
              im = V::ma(im, z.second, T(1));
            }
            V::store(y + size_t(q) * rhs + r, re, im);
          }
          for (; r < rhs; ++r) {
            Cx<T> sum{};
            for (uint32_t v = 0; v < nvec; ++v)
              sum = plus(sum, cm(p[size_t(q) * nvec + v],
                                 x[size_t(base + v) * rhs + r]));
            y[size_t(q) * rhs + r] =
                add ? plus(y[size_t(q) * rhs + r], sum) : sum;
          }
        }
      });
    }
  };

  std::vector<D> inverse_odd(const CloverSystem &c) {
    std::vector<D> out(c.oddA.size());
    for (size_t i = 0; i < c.oddA.size() / 36; ++i) {
      M6 a;
      for (int j = 0; j < 36; ++j) {
        auto z = c.oddA[i * 36 + j];
        a[j] = {z.re, z.im};
      }
      double pivot = 1e100;
      auto inv = invert_cholesky(a, pivot);
      for (int j = 0; j < 36; ++j)
        out[i * 36 + j] = {inv[j].real(), inv[j].imag()};
    }
    return out;
  }
  struct CPUMultigrid {
    const CloverSystem &c;
    const TwoLevel &h;
    CPUCloverSchur<float> &op;
    int workers;
    bool coarse_enabled = true;
    std::vector<Cx<float>> inv, residual, work, smooth, coarse_rhs, coarse_x;
    uint64_t calls = 0, schur_calls = 0;
    std::function<void(const Cx<float> *, Cx<float> *)> coarse_override;
    CPUMultigrid(const CloverSystem &c, const TwoLevel &h,
                 CPUCloverSchur<float> &op, int workers)
        : c(c), h(h), op(op), workers(workers),
          inv(block_cast<float>(inverse_odd(c))),
          residual(c.parity.oe.p.outputs()), work(residual.size()),
          smooth(residual.size()), coarse_rhs(size_t(h.coarse) * c.w.a.p.rhs),
          coarse_x(coarse_rhs.size()) {}
    void apply(const Cx<float> *r, Cx<float> *z) {
      ++calls;
      uint32_t rows = c.parity.oe.p.rows, rhs = c.w.a.p.rhs;
      float omega = float(h.omega);
      auto elements = [&](const auto &f) {
        chunks(rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j)
            f(j);
        });
      };
      auto refresh = [&]() {
        op.apply(z, work.data());
        ++schur_calls;
        elements([&](size_t j) { residual[j] = minus_z(r[j], work[j]); });
      };
      auto smooth_step = [&](const Cx<float> *input, bool first) {
        local_blocks(inv, input, first ? z : smooth.data(), rows, rhs, workers);
        elements([&](size_t j) {
          if (first) {
            z[j].re *= omega;
            z[j].im *= omega;
          } else {
            z[j].re += omega * smooth[j].re;
            z[j].im += omega * smooth[j].im;
          }
        });
      };
      if (!coarse_enabled) {
        for (int i = 0; i < h.pre_smooth + h.post_smooth; ++i) {
          if (i)
            refresh();
          smooth_step(i ? residual.data() : r, i == 0);
        }
        return;
      }
      for (int i = 0; i < h.pre_smooth; ++i) {
        if (i)
          refresh();
        smooth_step(i ? residual.data() : r, i == 0);
      }
      if (h.pre_smooth)
        refresh();
      h.restrict_vector(h.pre_smooth ? residual.data() : r, coarse_rhs.data(),
                        rhs);
      Cx<float> one{1, 0}, zero{};
      if (coarse_override)
        coarse_override(coarse_rhs.data(), coarse_x.data());
      else
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, int(h.coarse),
                    int(rhs), int(h.coarse),
                    blas_ptr(reinterpret_cast<Z *>(&one)),
                    blas_ptr(reinterpret_cast<const Z *>(h.inverse32.data())),
                    int(h.coarse),
                    blas_ptr(reinterpret_cast<const Z *>(coarse_rhs.data())),
                    int(rhs), blas_ptr(reinterpret_cast<Z *>(&zero)),
                    blas_ptr(reinterpret_cast<Z *>(coarse_x.data())), int(rhs));
      h.prolong(coarse_x.data(), z, rhs, h.pre_smooth > 0);
      for (int i = 0; i < h.post_smooth; ++i) {
        refresh();
        smooth_step(residual.data(), false);
      }
    }
  };

  struct GPUMultigrid {
    const CloverSystem &c;
    const TwoLevel &h;
    GPUCloverSchur &op;
    bool coarse_enabled = true;
    id<MTLBuffer> inv, basis, offsets, sites, aggregate, coarse_inverse,
        residual, work, smooth, coarse_rhs, coarse_x;
    std::array<id<MTLComputePipelineState>, 4> pipes;
    struct Params {
      uint32_t rows, coarse, nvec, mode;
      float omega;
    } params;
    uint64_t calls = 0, schur_calls = 0;
    std::function<void(id<MTLComputeCommandEncoder>, id<MTLBuffer>,
                       id<MTLBuffer>)>
        coarse_override;
    std::function<void()> coarse_check;
    void check_coarse() {
      if (coarse_check)
        coarse_check();
    }
    GPUMultigrid(const CloverSystem &c, const TwoLevel &h, GPUCloverSchur &op)
        : c(c), h(h), op(op) {
      auto put = [&](const void *data, size_t bytes) {
        auto b = [op.m.device newBufferWithLength:bytes
                                          options:MTLResourceStorageModeShared];
        if (!b)
          throw std::runtime_error("MG allocation failed");
        if (data)
          std::memcpy([b contents], data, bytes);
        return b;
      };
      auto inv32 = block_cast<float>(inverse_odd(c));
      inv = put(inv32.data(), inv32.size() * 8);
      basis = put(h.basis32.data(), h.basis32.size() * 8);
      offsets = put(h.offsets.data(), h.offsets.size() * 4);
      sites = put(h.sites.data(), h.sites.size() * 4);
      aggregate = put(h.aggregate.data(), h.aggregate.size() * 4);
      coarse_inverse = put(h.inverse32.empty() ? nullptr : h.inverse32.data(),
                           std::max(size_t(8), h.inverse32.size() * 8));
      size_t n = c.parity.oe.p.outputs() * 8;
      residual = put(nullptr, n);
      work = put(nullptr, n);
      smooth = put(nullptr, n);
      coarse_rhs = put(nullptr, size_t(h.coarse) * c.w.a.p.rhs * 8);
      coarse_x = put(nullptr, size_t(h.coarse) * c.w.a.p.rhs * 8);
      params = {h.rows, h.coarse, h.nvec, 0, float(h.omega)};
      auto values = [MTLFunctionConstantValues new];
      uint32_t rhs = c.w.a.p.rhs;
      [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
      const char *names[] = {"mg_scale_add", "mg_restrict", "mg_coarse",
                             "mg_prolong"};
      for (int j = 0; j < 4; ++j) {
        NSError *error = nil;
        auto f = [op.m.library
            newFunctionWithName:[NSString stringWithUTF8String:names[j]]
                 constantValues:values
                          error:&error];
        if (!f)
          throw std::runtime_error([[error localizedDescription] UTF8String]);
        pipes[j] =
            [op.m.device newComputePipelineStateWithFunction:f error:&error];
        if (!pipes[j])
          throw std::runtime_error([[error localizedDescription] UTF8String]);
      }
    }
    void dispatch(id<MTLComputeCommandEncoder> enc, int pipe,
                  std::initializer_list<id<MTLBuffer>> bs, size_t n,
                  uint32_t mode = 0) {
      params.mode = mode;
      [enc setComputePipelineState:pipes[pipe]];
      size_t j = 0;
      for (auto b : bs) {
        [enc setBuffer:b offset:0 atIndex:j];
        ++j;
      }
      [enc setBytes:&params length:sizeof(params) atIndex:j];
      [enc dispatchThreads:MTLSizeMake(n, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void encode(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> r,
                id<MTLBuffer> z) {
      ++calls;
      size_t n = c.parity.oe.p.outputs(), nc = size_t(h.coarse) * c.w.a.p.rhs;
      auto refresh = [&]() {
        op.encode(enc, z, work);
        ++schur_calls;
        op.difference(enc, r, work, residual);
      };
      auto smooth_step = [&](id<MTLBuffer> input, bool first) {
        op.block(enc, inv, input, first ? z : smooth);
        if (first)
          dispatch(enc, 0, {z, z}, n);
        else
          dispatch(enc, 0, {smooth, z}, n, 1);
      };
      if (!coarse_enabled) {
        for (int i = 0; i < h.pre_smooth + h.post_smooth; ++i) {
          if (i)
            refresh();
          smooth_step(i ? residual : r, i == 0);
        }
        return;
      }
      for (int i = 0; i < h.pre_smooth; ++i) {
        if (i)
          refresh();
        smooth_step(i ? residual : r, i == 0);
      }
      if (h.pre_smooth)
        refresh();
      dispatch(enc, 1,
               {h.pre_smooth ? residual : r, basis, offsets, sites, coarse_rhs},
               nc);
      if (coarse_override)
        coarse_override(enc, coarse_rhs, coarse_x);
      else
        dispatch(enc, 2, {coarse_inverse, coarse_rhs, coarse_x}, nc);
      dispatch(enc, 3, {coarse_x, basis, aggregate, z}, n,
               h.pre_smooth > 0 ? 1 : 0);
      for (int i = 0; i < h.post_smooth; ++i) {
        refresh();
        smooth_step(residual, false);
      }
    }
  };

} // namespace MetalCloverEnv
