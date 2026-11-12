#pragma once
#include "multigrid.hpp"

namespace MetalCloverEnv {

  struct ArnoldiData {
    std::vector<Cx<float>> h;
    std::vector<float> norm2;
  };
  struct SmallQR {
    uint32_t rhs, m;
    std::vector<CD> h, g, s;
    std::vector<double> c;
    std::vector<int> active, used;
    SmallQR(uint32_t rhs, uint32_t m)
        : rhs(rhs), m(m), h(size_t(rhs) * (m + 1) * m),
          g(size_t(rhs) * (m + 1)), s(size_t(rhs) * m), c(size_t(rhs) * m),
          active(rhs), used(rhs) {}
    CD &H(uint32_t r, uint32_t i, uint32_t j) {
      return h[size_t(r) * (m + 1) * m + i * m + j];
    }
    void begin(const std::vector<float> &current,
               const std::vector<float> &initial, double tol) {
      std::fill(h.begin(), h.end(), CD{});
      std::fill(g.begin(), g.end(), CD{});
      std::fill(used.begin(), used.end(), 0);
      for (uint32_t r = 0; r < rhs; ++r) {
        active[r] = current[r] > tol * tol * initial[r];
        g[size_t(r) * (m + 1)] = std::sqrt(current[r]);
      }
    }
    std::vector<float> factors(const std::vector<float> &norm) {
      std::vector<float> f(rhs);
      for (uint32_t r = 0; r < rhs; ++r) {
        if (active[r]) {
          if (!(norm[r] > 0) || !std::isfinite(norm[r]))
            throw std::runtime_error("FGMRES Arnoldi breakdown");
          f[r] = 1 / std::sqrt(norm[r]);
        }
      }
      return f;
    }
    void update(uint32_t k, const ArnoldiData &a,
                const std::vector<float> &initial, double tol) {
      for (uint32_t r = 0; r < rhs; ++r)
        if (active[r]) {
          for (uint32_t j = 0; j <= k; ++j) {
            auto z = a.h[size_t(j) * rhs + r];
            H(r, j, k) = {z.re, z.im};
          }
          H(r, k + 1, k) = std::sqrt(a.norm2[r]);
          for (uint32_t j = 0; j < k; ++j) {
            CD x = H(r, j, k), y = H(r, j + 1, k), sn = s[size_t(r) * m + j];
            double cs = c[size_t(r) * m + j];
            H(r, j, k) = cs * x + sn * y;
            H(r, j + 1, k) = -std::conj(sn) * x + cs * y;
          }
          CD f = H(r, k, k), b = H(r, k + 1, k);
          double n = std::hypot(std::abs(f), std::abs(b));
          if (!(n > 0) || !std::isfinite(n))
            throw std::runtime_error("FGMRES small QR breakdown");
          CD alpha = std::abs(f) > 0 ? f / std::abs(f) : CD(1);
          double cs = std::abs(f) / n;
          CD sn = alpha * std::conj(b) / n;
          c[size_t(r) * m + k] = cs;
          s[size_t(r) * m + k] = sn;
          H(r, k, k) = alpha * n;
          H(r, k + 1, k) = 0;
          CD x = g[size_t(r) * (m + 1) + k], y = g[size_t(r) * (m + 1) + k + 1];
          g[size_t(r) * (m + 1) + k] = cs * x + sn * y;
          g[size_t(r) * (m + 1) + k + 1] = -std::conj(sn) * x + cs * y;
          used[r] = int(k + 1);
          if (std::norm(g[size_t(r) * (m + 1) + k + 1]) <=
              tol * tol * initial[r])
            active[r] = 0;
        }
    }
    std::vector<Cx<float>> solution(uint32_t k) {
      std::vector<Cx<float>> y(size_t(k) * rhs);
      for (uint32_t r = 0; r < rhs; ++r) {
        std::vector<CD> v(used[r]);
        for (int i = used[r] - 1; i >= 0; --i) {
          CD z = g[size_t(r) * (m + 1) + i];
          for (int j = i + 1; j < used[r]; ++j)
            z -= H(r, uint32_t(i), uint32_t(j)) * v[j];
          v[i] = z / H(r, uint32_t(i), uint32_t(i));
          if (!std::isfinite(v[i].real()) || !std::isfinite(v[i].imag()))
            throw std::runtime_error("nonfinite FGMRES solution");
        }
        for (int i = 0; i < used[r]; ++i)
          y[size_t(i) * rhs + r] = {float(v[i].real()), float(v[i].imag())};
      }
      return y;
    }
  };

  template <class Backend>
  InnerInfo fgmres_iterations(Backend &b, double tol, int limit) {
    auto initial = b.current_norms(), current = initial;
    SmallQR qr(b.rhs, b.restart);
    std::vector<int> counts(b.rhs);
    int total = 0;
    auto done = [&]() {
      for (uint32_t r = 0; r < b.rhs; ++r) {
        if (!std::isfinite(current[r]))
          throw std::runtime_error("FGMRES nonfinite true FP32 residual");
        if (current[r] > tol * tol * initial[r])
          return false;
      }
      return true;
    };
    while (!done()) {
      if (total >= limit) {
        if (b.allow_inexact)
          break;
        throw std::runtime_error("FGMRES iteration limit");
      }
      qr.begin(current, initial, tol);
      b.normalize(0, qr.factors(current), true);
      uint32_t k = 0;
      for (uint32_t j = 0; j < b.restart && total < limit; ++j) {
        for (uint32_t r = 0; r < b.rhs; ++r)
          if (qr.active[r])
            ++counts[r];
        auto a = b.step(j);
        qr.update(j, a, initial, tol);
        ++total;
        k = j + 1;
        if (std::none_of(qr.active.begin(), qr.active.end(),
                         [](int v) { return v != 0; }))
          break;
        if (j + 1 < b.restart)
          b.normalize(j + 1, qr.factors(a.norm2), false);
      }
      b.combine(k, qr.solution(k));
      current = b.refresh();
      if (b.progress && total % 64 == 0) {
        double worst = 0;
        for (uint32_t r = 0; r < b.rhs; ++r)
          if (initial[r] > 0)
            worst = std::max(worst, std::sqrt(double(current[r]) / initial[r]));
        b.progress(total, worst);
      }
    }
    return {*std::min_element(counts.begin(), counts.end()),
            *std::max_element(counts.begin(), counts.end()), total};
  }

  template <class Operator, class Preconditioner> struct CPUGenericFGMRES {
    Operator &op;
    Preconditioner &pre;
    int workers;
    uint32_t rhs, restart, components;
    size_t n;
    std::vector<Cx<float>> saved, r, w, delta, vcur, zcur, V, Z, partial;
    uint64_t fine_calls = 0;
    bool allow_inexact = false;
    std::function<void(int, double)> progress;
    CPUGenericFGMRES(uint32_t components, uint32_t rhs, Operator &op,
                     Preconditioner &pre, int workers, uint32_t restart)
        : op(op), pre(pre), workers(workers), rhs(rhs), restart(restart),
          components(components), n(size_t(components) * rhs), saved(n), r(n),
          w(n), delta(n), vcur(n), zcur(n), V((restart + 1) * n),
          Z(restart * n), partial(size_t(workers) * restart * rhs) {}
    CPUGenericFGMRES(const CloverSystem &c, Operator &op, Preconditioner &pre,
                     int workers, uint32_t restart)
        : CPUGenericFGMRES(c.parity.oe.p.rows * 12, c.w.a.p.rhs, op, pre,
                           workers, restart) {}
    template <class F> void elements(const F &f) {
      chunks(components, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t q = first; q < last; ++q)
          for (uint32_t k = 0; k < rhs; ++k)
            f(q * rhs + k, k);
      });
    }
    std::vector<Cx<float>> dots(uint32_t k) {
      chunks(
          components, workers, [&](int worker, uint32_t first, uint32_t last) {
            auto *sums = partial.data() + size_t(worker) * restart * rhs;
            std::fill(sums, sums + size_t(k) * rhs, Cx<float>{});
            for (uint32_t j = 0; j < k; ++j)
              for (size_t q = first; q < last; ++q)
                for (uint32_t v = 0; v < rhs; ++v) {
                  auto a = V[size_t(j) * n + q * rhs + v], b = w[q * rhs + v];
                  auto &z = sums[size_t(j) * rhs + v];
                  z.re += a.re * b.re + a.im * b.im;
                  z.im += a.re * b.im - a.im * b.re;
                }
          });
      std::vector<Cx<float>> out(size_t(k) * rhs);
      for (int worker = 0; worker < workers; ++worker)
        for (size_t j = 0; j < out.size(); ++j)
          out[j] = plus(out[j], partial[size_t(worker) * restart * rhs + j]);
      return out;
    }
    std::vector<float> norm(const std::vector<Cx<float>> &x) {
      chunks(components, workers,
             [&](int worker, uint32_t first, uint32_t last) {
               auto *sums = partial.data() + size_t(worker) * restart * rhs;
               std::fill(sums, sums + rhs, Cx<float>{});
               for (size_t q = first; q < last; ++q)
                 for (uint32_t v = 0; v < rhs; ++v) {
                   auto z = x[q * rhs + v];
                   sums[v].re += z.re * z.re + z.im * z.im;
                 }
             });
      std::vector<float> out(rhs);
      for (int worker = 0; worker < workers; ++worker)
        for (uint32_t v = 0; v < rhs; ++v)
          out[v] += partial[size_t(worker) * restart * rhs + v].re;
      return out;
    }
    std::vector<float> current_norms() { return norm(r); }
    void normalize(uint32_t k, const std::vector<float> &f, bool initial) {
      const auto &x = initial ? r : w;
      elements([&](size_t j, uint32_t v) {
        V[size_t(k) * n + j] = {x[j].re * f[v], x[j].im * f[v]};
      });
    }
    ArnoldiData step(uint32_t k) {
      elements([&](size_t j, uint32_t) { vcur[j] = V[size_t(k) * n + j]; });
      pre.apply(vcur.data(), zcur.data());
      elements([&](size_t j, uint32_t) { Z[size_t(k) * n + j] = zcur[j]; });
      op.apply(zcur.data(), w.data());
      ++fine_calls;
      std::vector<Cx<float>> h;
      for (int pass = 0; pass < 2; ++pass) {
        auto a = dots(k + 1);
        if (pass == 0)
          h = a;
        else
          for (size_t j = 0; j < h.size(); ++j)
            h[j] = plus(h[j], a[j]);
        elements([&](size_t j, uint32_t v) {
          auto sum = w[j];
          for (uint32_t i = 0; i <= k; ++i)
            sum =
                minus_z(sum, cm(V[size_t(i) * n + j], a[size_t(i) * rhs + v]));
          w[j] = sum;
        });
      }
      return {h, norm(w)};
    }
    void combine(uint32_t k, const std::vector<Cx<float>> &y) {
      elements([&](size_t j, uint32_t v) {
        auto z = delta[j];
        for (uint32_t i = 0; i < k; ++i)
          z = plus(z, cm(Z[size_t(i) * n + j], y[size_t(i) * rhs + v]));
        delta[j] = z;
      });
    }
    std::vector<float> refresh() {
      op.apply(delta.data(), w.data());
      ++fine_calls;
      elements([&](size_t j, uint32_t) { r[j] = minus_z(saved[j], w[j]); });
      return norm(r);
    }
    template <class T>
    InnerInfo run(const std::vector<Cx<T>> &input, double tol, int limit) {
      elements([&](size_t j, uint32_t) {
        r[j] = saved[j] = {float(input[j].re), float(input[j].im)};
        delta[j] = {};
      });
      return fgmres_iterations(*this, tol, limit);
    }
    const Cx<float> *correction() const { return delta.data(); }
  };

  using CPUFGMRES = CPUGenericFGMRES<CPUCloverSchur<float>, CPUMultigrid>;

  template <class Operator, class Preconditioner> struct GPUGenericFGMRES {
    Operator &op;
    Preconditioner &pre;
    uint32_t rhs, restart;
    size_t n;
    id<MTLBuffer> saved, r, w, delta, vcur, zcur, V, Z, partial, coefficients,
        h, partial_norm, norms, factors, y;
    std::array<id<MTLComputePipelineState>, 12> pipes;
    struct Params {
      uint32_t components, parts, k, stage;
    } params;
    uint64_t fine_calls = 0;
    bool allow_inexact = false;
    std::function<void(int, double)> progress;
    std::function<void(id<MTLBuffer>, id<MTLBuffer>)> blocking_pre;
    GPUGenericFGMRES(uint32_t components, uint32_t rhs, Operator &op,
                     Preconditioner &pre, uint32_t restart)
        : op(op), pre(pre), rhs(rhs), restart(restart),
          n(size_t(components) * rhs) {
      auto alloc = [&](size_t bytes) {
        auto b = [op.m.device newBufferWithLength:bytes
                                          options:MTLResourceStorageModeShared];
        if (!b)
          throw std::runtime_error("FGMRES allocation failed");
        return b;
      };
      saved = alloc(n * 8);
      r = alloc(n * 8);
      w = alloc(n * 8);
      delta = alloc(n * 8);
      vcur = alloc(n * 8);
      zcur = alloc(n * 8);
      V = alloc((restart + 1) * n * 8);
      Z = alloc(restart * n * 8);
      params = {components, (components + 255) / 256, 0, 0};
      partial = alloc(size_t(restart) * params.parts * rhs * 8);
      coefficients = alloc(restart * rhs * 8);
      h = alloc(restart * rhs * 8);
      partial_norm = alloc(size_t(params.parts) * rhs * 4);
      norms = alloc(rhs * 4);
      factors = alloc(rhs * 4);
      y = alloc(restart * rhs * 8);
      const char *names[] = {
          "gm_init",        "gm_store",   "gm_extract",       "gm_normalize",
          "gm_dots",        "gm_reduce",  "gm_orthogonalize", "gm_norm_partial",
          "gm_norm_reduce", "gm_combine", "clover_subtract",  "mg_scale_add"};
      auto values = [MTLFunctionConstantValues new];
      [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
      for (int j = 0; j < 12; ++j) {
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
    GPUGenericFGMRES(const CloverSystem &c, Operator &op, Preconditioner &pre,
                     uint32_t restart)
        : GPUGenericFGMRES(c.parity.oe.p.rows * 12, c.w.a.p.rhs, op, pre,
                           restart) {}
    void dispatch(id<MTLComputeCommandEncoder> enc, int p,
                  std::initializer_list<id<MTLBuffer>> bs, size_t size) {
      [enc setComputePipelineState:pipes[p]];
      size_t j = 0;
      for (auto b : bs) {
        [enc setBuffer:b offset:0 atIndex:j];
        ++j;
      }
      [enc setBytes:&params length:sizeof(params) atIndex:j];
      [enc dispatchThreads:MTLSizeMake(size, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void finish(id<MTLCommandBuffer> cmd, id<MTLComputeCommandEncoder> enc) {
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
      if (cmd.status == MTLCommandBufferStatusError)
        throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);
    }
    void encode_norm(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x) {
      dispatch(enc, 7, {x, partial_norm}, size_t(params.parts) * rhs);
      dispatch(enc, 8, {partial_norm, norms}, rhs);
    }
    std::vector<float> read_norm() {
      auto *p = static_cast<float *>([norms contents]);
      return {p, p + rhs};
    }
    std::vector<float> current_norms() {
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        encode_norm(enc, r);
        finish(cmd, enc);
        return read_norm();
      }
    }
    void normalize(uint32_t k, const std::vector<float> &f, bool initial) {
      std::memcpy([factors contents], f.data(), rhs * 4);
      params.k = k;
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        dispatch(enc, 3, {initial ? r : w, V, factors}, n);
        finish(cmd, enc);
      }
    }
    ArnoldiData step(uint32_t k) {
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        params.k = k;
        dispatch(enc, 2, {V, vcur}, n);
        if (blocking_pre) {
          finish(cmd, enc);
          blocking_pre(vcur, zcur);
          cmd = [op.m.queue commandBuffer];
          enc = [cmd computeCommandEncoder];
        } else
          pre.encode(enc, vcur, zcur);
        params.k = k;
        params.stage = 0;
        dispatch(enc, 1, {zcur, Z}, n);
        op.encode(enc, zcur, w);
        ++fine_calls;
        params.k = k + 1;
        for (uint32_t pass = 0; pass < 2; ++pass) {
          params.stage = pass;
          dispatch(enc, 4, {V, w, partial}, size_t(k + 1) * params.parts * rhs);
          dispatch(enc, 5, {partial, coefficients, h}, size_t(k + 1) * rhs);
          dispatch(enc, 6, {V, coefficients, w}, n);
        }
        encode_norm(enc, w);
        finish(cmd, enc);
        pre.check_coarse();
        auto *hp = static_cast<Cx<float> *>([h contents]);
        return {{hp, hp + size_t(k + 1) * rhs}, read_norm()};
      }
    }
    void combine(uint32_t k, const std::vector<Cx<float>> &values) {
      std::memcpy([y contents], values.data(), values.size() * 8);
      params.k = k;
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        dispatch(enc, 9, {Z, y, delta}, n);
        finish(cmd, enc);
      }
    }
    std::vector<float> refresh() {
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        op.encode(enc, delta, w);
        ++fine_calls;
        op.difference(enc, saved, w, r);
        encode_norm(enc, r);
        finish(cmd, enc);
        return read_norm();
      }
    }
    InnerInfo run_device(id<MTLBuffer> input, double tol, int limit) {
      @autoreleasepool {
        auto cmd = [op.m.queue commandBuffer];
        auto enc = [cmd computeCommandEncoder];
        dispatch(enc, 0, {input, saved, r, delta}, n);
        finish(cmd, enc);
      }
      return fgmres_iterations(*this, tol, limit);
    }
    const Cx<float> *correction() const {
      return static_cast<Cx<float> *>([delta contents]);
    }
  };
  using GPUFGMRES = GPUGenericFGMRES<GPUCloverSchur, GPUMultigrid>;

} // namespace MetalCloverEnv
