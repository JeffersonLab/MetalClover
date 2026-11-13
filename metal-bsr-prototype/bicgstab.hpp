#pragma once
#include "solve_common.hpp"

namespace MetalCloverEnv {

  template <class T> Cx<T> minus_z(Cx<T> a, Cx<T> b) {
    return {a.re - b.re, a.im - b.im};
  }
  template <class T> Cx<T> divide_z(Cx<T> a, Cx<T> b) {
    if (std::abs(b.re) >= std::abs(b.im)) {
      T q = b.im / b.re, d = b.re + b.im * q;
      return {(a.re + a.im * q) / d, (a.im - a.re * q) / d};
    }
    T q = b.re / b.im, d = b.im + b.re * q;
    return {(a.re * q + a.im) / d, (a.im * q - a.re) / d};
  }
  template <class T> bool valid_z(Cx<T> z) {
    return std::isfinite(z.re) && std::isfinite(z.im) &&
           std::max(std::abs(z.re), std::abs(z.im)) > 0;
  }

  template <class T> struct CPUBiCGStab {
    const System &a;
    int workers;
    uint32_t components;
    std::vector<Cx<T>> delta, r, shadow, p, v, t, partial, rho, old, alpha,
        omega, beta, z, w, initial;
    std::vector<int> active, iterations;
    std::function<void(const Cx<T> *, Cx<T> *)> apply_override;
    bool allow_inexact = false;
    bool stable_reductions = false;
    std::function<void(int, double)> progress;
    CPUBiCGStab(const System &a, int workers)
        : a(a), workers(workers), components(a.p.rows * a.p.colors * a.p.spins),
          delta(a.p.outputs()), r(delta.size()), shadow(delta.size()),
          p(delta.size()), v(delta.size()), t(delta.size()),
          partial(workers * a.p.rhs), rho(a.p.rhs), old(rho), alpha(rho),
          omega(rho), beta(rho), z(rho), w(rho), initial(rho), active(a.p.rhs),
          iterations(a.p.rhs) {}
    void dot(const std::vector<Cx<T>> &x, const std::vector<Cx<T>> &y,
             std::vector<Cx<T>> &out) {
      uint32_t rhs = a.p.rhs;
      if (stable_reductions) {
        uint32_t parts = (components + 255) / 256;
        partial.resize(size_t(parts) * rhs);
        chunks(parts, workers, [&](int, uint32_t first, uint32_t last) {
          for (uint32_t part = first; part < last; ++part) {
            auto *sum = partial.data() + size_t(part) * rhs;
            std::fill(sum, sum + rhs, Cx<T>{});
            for (uint32_t q = part * 256;
                 q < std::min((part + 1) * 256, components); ++q)
              for (uint32_t r = 0; r < rhs; ++r) {
                auto xx = x[size_t(q) * rhs + r], yy = y[size_t(q) * rhs + r];
                sum[r].re = std::fma(xx.re, yy.re, sum[r].re);
                sum[r].re = std::fma(xx.im, yy.im, sum[r].re);
                sum[r].im = std::fma(xx.re, yy.im, sum[r].im);
                sum[r].im = std::fma(-xx.im, yy.re, sum[r].im);
              }
          }
        });
        for (uint32_t r = 0; r < rhs; ++r) {
          double re = 0, im = 0;
          for (uint32_t p = 0; p < parts; ++p) {
            auto z = partial[size_t(p) * rhs + r];
            re += z.re;
            im += z.im;
          }
          out[r] = {T(re), T(im)};
        }
        return;
      }
      chunks(components, workers,
             [&](int worker, uint32_t first, uint32_t last) {
               auto *sums = partial.data() + worker * rhs;
               std::fill(sums, sums + rhs, Cx<T>{});
               for (size_t j = first; j < last; ++j)
                 for (uint32_t k = 0; k < rhs; ++k) {
                   auto xx = x[j * rhs + k], yy = y[j * rhs + k];
                   sums[k].re += xx.re * yy.re + xx.im * yy.im;
                   sums[k].im += xx.re * yy.im - xx.im * yy.re;
                 }
             });
      std::fill(out.begin(), out.end(), Cx<T>{});
      for (int j = 0; j < workers; ++j)
        for (uint32_t k = 0; k < rhs; ++k)
          out[k] = plus(out[k], partial[j * rhs + k]);
    }
    template <class F> void vectors(const F &f) {
      uint32_t rhs = a.p.rhs;
      chunks(components, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t j = first; j < last; ++j)
          for (uint32_t k = 0; k < rhs; ++k)
            f(j * rhs + k, k);
      });
    }
    void apply(const Cx<T> *x, Cx<T> *y) {
      if (apply_override)
        apply_override(x, y);
      else
        apply_cpu(a, x, y, workers);
    }
    template <class R>
    InnerInfo run(const std::vector<Cx<R>> &residual, double tol, int limit) {
      uint32_t rhs = a.p.rhs;
      vectors([&](size_t j, uint32_t) {
        r[j] = shadow[j] = {T(residual[j].re), T(residual[j].im)};
        delta[j] = p[j] = v[j] = {};
      });
      dot(r, r, initial);
      for (uint32_t k = 0; k < rhs; ++k) {
        old[k] = alpha[k] = omega[k] = {1, 0};
        active[k] = initial[k].re > 0;
        iterations[k] = 0;
      }
      int it = 0;
      while (std::any_of(active.begin(), active.end(),
                         [](int k) { return k != 0; })) {
        if (it >= limit) {
          if (allow_inexact)
            break;
          throw std::runtime_error("CPU BiCGStab iteration limit");
        }
        dot(shadow, r, rho);
        for (uint32_t k = 0; k < rhs; ++k)
          if (active[k]) {
            if (!valid_z(rho[k]) || !valid_z(old[k]) || !valid_z(omega[k]))
              throw std::runtime_error("CPU BiCGStab rho/omega breakdown");
            beta[k] =
                cm(divide_z(rho[k], old[k]), divide_z(alpha[k], omega[k]));
          }
        vectors([&](size_t j, uint32_t k) {
          p[j] =
              active[k]
                  ? plus(r[j], cm(beta[k], minus_z(p[j], cm(omega[k], v[j]))))
                  : Cx<T>{};
        });
        apply(p.data(), v.data());
        dot(shadow, v, z);
        for (uint32_t k = 0; k < rhs; ++k)
          if (active[k]) {
            if (!valid_z(z[k]))
              throw std::runtime_error("CPU BiCGStab alpha breakdown");
            alpha[k] = divide_z(rho[k], z[k]);
          }
        vectors([&](size_t j, uint32_t k) {
          if (active[k]) {
            delta[j] = plus(delta[j], cm(alpha[k], p[j]));
            r[j] = minus_z(r[j], cm(alpha[k], v[j]));
          }
        });
        dot(r, r, z);
        for (uint32_t k = 0; k < rhs; ++k)
          if (active[k] && z[k].re <= T(tol * tol) * initial[k].re) {
            active[k] = 0;
            ++iterations[k];
          }
        if (std::any_of(active.begin(), active.end(),
                        [](int k) { return k != 0; })) {
          apply(r.data(), t.data());
          dot(t, r, z);
          dot(t, t, w);
          for (uint32_t k = 0; k < rhs; ++k)
            if (active[k]) {
              if (!(w[k].re > 0) || !std::isfinite(w[k].re))
                throw std::runtime_error("CPU BiCGStab t breakdown");
              omega[k] = {z[k].re / w[k].re, z[k].im / w[k].re};
              if (!valid_z(omega[k]))
                throw std::runtime_error("CPU BiCGStab zero omega");
            }
          vectors([&](size_t j, uint32_t k) {
            if (active[k]) {
              delta[j] = plus(delta[j], cm(omega[k], r[j]));
              r[j] = minus_z(r[j], cm(omega[k], t[j]));
            }
          });
          dot(r, r, z);
          for (uint32_t k = 0; k < rhs; ++k)
            if (active[k]) {
              ++iterations[k];
              old[k] = rho[k];
              if (!std::isfinite(z[k].re))
                throw std::runtime_error("CPU BiCGStab nonfinite norm");
              if (z[k].re <= T(tol * tol) * initial[k].re)
                active[k] = 0;
            }
        }
        ++it;
        if (progress && it % 64 == 0) {
          double worst = 0;
          for (uint32_t r = 0; r < rhs; ++r)
            if (initial[r].re > 0)
              worst =
                  std::max(worst, std::sqrt(double(z[r].re) / initial[r].re));
          progress(it, worst);
        }
      }
      return {*std::min_element(iterations.begin(), iterations.end()),
              *std::max_element(iterations.begin(), iterations.end()), it};
    }
    const Cx<T> *correction() const { return delta.data(); }
  };

  struct GPUBiCGStab {
    GPUProblem matrix;
    MetalRun &m;
    std::array<id<MTLComputePipelineState>, 5> pipes;
    id<MTLBuffer> r, shadow, t, delta, partial, second, scalars;
    struct Params {
      uint32_t components, parts, limit;
      float tolerance2;
      uint32_t stage;
    } params;
    int chunk, workers;
    bool allow_inexact = false;
    uint32_t rhs;
    std::function<void(id<MTLComputeCommandEncoder>, id<MTLBuffer>,
                       id<MTLBuffer>)>
        encode_override;
    id<MTLBuffer> device_initial = nil;
    std::function<void(int, double)> progress;
    GPUBiCGStab(MetalRun &m, const System &a, int chunk, int workers)
        : matrix(m, a.p), m(m), chunk(chunk), workers(workers), rhs(a.p.rhs) {
      params = {a.p.rows * 12, (a.p.rows * 12 + 255) / 256, 0, 0, 0};
      auto values = [MTLFunctionConstantValues new];
      [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
      const char *names[] = {"bi_partial", "bi_scalar", "bi_direction",
                             "bi_update", "bi_initial"};
      NSError *error = nil;
      for (int j = 0; j < 5; ++j) {
        auto fn = [m.library
            newFunctionWithName:[NSString stringWithUTF8String:names[j]]
                 constantValues:values
                          error:&error];
        if (!fn)
          throw std::runtime_error([[error localizedDescription] UTF8String]);
        pipes[j] =
            [m.device newComputePipelineStateWithFunction:fn error:&error];
        if (!pipes[j])
          throw std::runtime_error([[error localizedDescription] UTF8String]);
      }
      auto alloc = [&](size_t bytes) {
        auto b = [m.device newBufferWithLength:bytes
                                       options:MTLResourceStorageModeShared];
        if (!b)
          throw std::runtime_error("BiCGStab allocation failed");
        return b;
      };
      // The multiply's optional input may be empty; BiCGStab needs a full
      // direction workspace.
      if ([matrix.buffers[3] length] < a.p.outputs() * 8)
        matrix.buffers[3] = alloc(a.p.outputs() * 8);
      r = alloc(a.p.outputs() * 8);
      shadow = alloc(a.p.outputs() * 8);
      t = alloc(a.p.outputs() * 8);
      delta = alloc(a.p.outputs() * 8);
      partial = alloc(size_t(params.parts) * rhs * 8);
      second = alloc(size_t(params.parts) * rhs * 8);
      scalars = alloc(rhs * 10 * 8);
    }
    InnerInfo run(const std::vector<D> &residual, double tol, int limit) {
      device_initial = nil;
      auto *rr = static_cast<Cx<float> *>([r contents]);
      auto *sh = static_cast<Cx<float> *>([shadow contents]);
      auto *pp = static_cast<Cx<float> *>([matrix.buffers[3] contents]);
      auto *vv = static_cast<Cx<float> *>([matrix.buffers[4] contents]);
      auto *dd = static_cast<Cx<float> *>([delta contents]);
      chunks(
          params.components, workers, [&](int, uint32_t first, uint32_t last) {
            for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
              rr[j] = sh[j] = {float(residual[j].re), float(residual[j].im)};
              pp[j] = vv[j] = dd[j] = {};
            }
          });
      return iterate(tol, limit);
    }
    InnerInfo run_device(id<MTLBuffer> initial, double tol, int limit) {
      device_initial = initial;
      return iterate(tol, limit);
    }
    InnerInfo iterate(double tol, int limit) {
      params.limit = limit;
      params.tolerance2 = float(tol * tol);
      auto *ss = static_cast<Cx<float> *>([scalars contents]);
      int dispatched = 0;
      for (int offset = 0; offset < limit; offset += chunk) {
        @autoreleasepool {
          auto cmd = [m.queue commandBuffer];
          auto enc = [cmd computeCommandEncoder];
          auto dispatch = [&](int pipe, std::initializer_list<id<MTLBuffer>> bs,
                              size_t n, uint32_t stage = 0) {
            params.stage = stage;
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
          };
          auto dot = [&](id<MTLBuffer> x, id<MTLBuffer> y, id<MTLBuffer> out) {
            dispatch(0, {x, y, out}, size_t(params.parts) * rhs);
          };
          auto scalar = [&](uint32_t stage) {
            dispatch(1, {partial, second, scalars}, rhs, stage);
          };
          auto apply = [&](id<MTLBuffer> x, id<MTLBuffer> y) {
            if (encode_override) {
              encode_override(enc, x, y);
              return;
            }
            [enc setComputePipelineState:matrix.pipeline];
            for (size_t j = 0; j < matrix.buffers.size(); ++j)
              [enc setBuffer:matrix.buffers[j] offset:0 atIndex:j];
            [enc setBuffer:x offset:0 atIndex:3];
            [enc setBuffer:y offset:0 atIndex:4];
            [enc setBytes:&matrix.rows length:4 atIndex:7];
            [enc dispatchThreads:MTLSizeMake(matrix.outputs, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
            [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
          };
          if (offset == 0) {
            if (device_initial) {
              dispatch(4,
                       {device_initial, r, shadow, matrix.buffers[3],
                        matrix.buffers[4], delta},
                       matrix.outputs);
              device_initial = nil;
            }
            dot(r, r, partial);
            scalar(0);
          }
          for (int it = offset; it < std::min(limit, offset + chunk); ++it) {
            dot(shadow, r, partial);
            scalar(1);
            dispatch(2, {r, matrix.buffers[3], matrix.buffers[4], scalars},
                     matrix.outputs);
            apply(matrix.buffers[3], matrix.buffers[4]);
            dot(shadow, matrix.buffers[4], partial);
            scalar(2);
            dispatch(
                3, {delta, r, matrix.buffers[3], matrix.buffers[4], t, scalars},
                matrix.outputs, 0);
            dot(r, r, partial);
            scalar(3);
            apply(r, t);
            dot(t, r, partial);
            dot(t, t, second);
            scalar(4);
            dispatch(
                3, {delta, r, matrix.buffers[3], matrix.buffers[4], t, scalars},
                matrix.outputs, 1);
            dot(r, r, partial);
            scalar(5);
            ++dispatched;
          }
          [enc endEncoding];
          [cmd commit];
          [cmd waitUntilCompleted];
          if (cmd.status == MTLCommandBufferStatusError)
            throw std::runtime_error(
                [[cmd.error localizedDescription] UTF8String]);
        }
        bool active = false;
        for (uint32_t k = 0; k < rhs; ++k) {
          if (ss[8 * rhs + k].re)
            throw std::runtime_error("GPU BiCGStab scalar breakdown");
          active |= ss[6 * rhs + k].re > 0;
        }
        if (progress && dispatched % 64 == 0) {
          double worst = 0;
          for (uint32_t r = 0; r < rhs; ++r)
            if (ss[5 * rhs + r].re > 0)
              worst = std::max(worst, std::sqrt(double(ss[9 * rhs + r].re) /
                                                ss[5 * rhs + r].re));
          progress(dispatched, worst);
        }
        if (!active)
          break;
      }
      int lo = limit, hi = 0;
      std::vector<int> counts(rhs);
      for (uint32_t k = 0; k < rhs; ++k) {
        if (!allow_inexact &&
            ss[9 * rhs + k].re > params.tolerance2 * ss[5 * rhs + k].re)
          throw std::runtime_error("GPU BiCGStab iteration limit");
        counts[k] = int(ss[7 * rhs + k].re);
        lo = std::min(lo, counts[k]);
        hi = std::max(hi, int(ss[7 * rhs + k].re));
      }
      return {lo, hi, dispatched, counts};
    }
    const Cx<float> *correction() const {
      return static_cast<Cx<float> *>([delta contents]);
    }
  };

} // namespace MetalCloverEnv
