#pragma once
#include "bsr_common.hpp"
#include <array>
#include <memory>

namespace MetalCloverEnv {

  template <class T> struct Cx {
    T re = 0, im = 0;
  };
  using D = Cx<double>;
  template <class T> Cx<T> cm(Cx<T> a, Cx<T> b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
  }
  template <class T> Cx<T> plus(Cx<T> a, Cx<T> b) {
    return {a.re + b.re, a.im + b.im};
  }
  template <class F> void chunks(uint32_t n, int workers, const F &f) {
    if (workers == 1) {
      f(0, 0, n);
      return;
    }
    dispatch_apply(workers,
                   dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t w) {
                     f(int(w), uint32_t(uint64_t(n) * w / workers),
                       uint32_t(uint64_t(n) * (w + 1) / workers));
                   });
  }

  template <class T> struct SIMD;
  template <> struct SIMD<float> {
    using V = float32x4_t;
    static constexpr int width = 4;
    static V zero() { return vdupq_n_f32(0); }
    static std::pair<V, V> load(const Cx<float> *p) {
      auto z = vld2q_f32(reinterpret_cast<const float *>(p));
      return {z.val[0], z.val[1]};
    }
    static V ma(V a, V x, float v) { return vfmaq_n_f32(a, x, v); }
    static V ms(V a, V x, float v) { return vfmsq_n_f32(a, x, v); }
    static void store(Cx<float> *p, V re, V im) {
      float32x4x2_t z = {{re, im}};
      vst2q_f32(reinterpret_cast<float *>(p), z);
    }
  };
  template <> struct SIMD<double> {
    using V = float64x2_t;
    static constexpr int width = 2;
    static V zero() { return vdupq_n_f64(0); }
    static std::pair<V, V> load(const Cx<double> *p) {
      auto z = vld2q_f64(reinterpret_cast<const double *>(p));
      return {z.val[0], z.val[1]};
    }
    static V ma(V a, V x, double v) { return vfmaq_n_f64(a, x, v); }
    static V ms(V a, V x, double v) { return vfmsq_n_f64(a, x, v); }
    static void store(Cx<double> *p, V re, V im) {
      float64x2x2_t z = {{re, im}};
      vst2q_f64(reinterpret_cast<double *>(p), z);
    }
  };

  struct System {
    Problem p;
    std::vector<Cx<float>> af;
    std::vector<D> ad, b, truth;
    std::vector<double> bnorm;
    double mass;
    // Optional on-site clover correction: two row-major 6x6 chiral blocks/site.
    std::vector<Cx<float>> local32;
    std::vector<D> local64;
    double solution_error_limit = 0;
    std::function<void(const Cx<float> *, Cx<float> *, int)> apply32;
    std::function<void(const D *, D *, int)> apply64;
  };

  template <class T, class Coefficient = T>
  void local_blocks(const std::vector<Cx<Coefficient>> &blocks, const Cx<T> *x,
                    Cx<T> *y, uint32_t rows, uint32_t rhs, int workers,
                    T beta = T(0)) {
    using V = SIMD<T>;
    chunks(rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (uint32_t i = first; i < last; ++i)
        for (uint32_t ch = 0; ch < 2; ++ch)
          for (uint32_t q = 0; q < 6; ++q) {
            size_t out = ((size_t(i) * 3 + q % 3) * 4 + ch * 2 + q / 3) * rhs;
            uint32_t r = 0;
            for (; r + V::width <= rhs; r += V::width) {
              auto re = V::zero(), im = V::zero();
              for (uint32_t k = 0; k < 6; ++k) {
                auto a = blocks[(size_t(i) * 2 + ch) * 36 + q * 6 + k];
                auto v = V::load(
                    x + ((size_t(i) * 3 + k % 3) * 4 + ch * 2 + k / 3) * rhs +
                    r);
                re = V::ma(re, v.first, a.re);
                re = V::ms(re, v.second, a.im);
                im = V::ma(im, v.second, a.re);
                im = V::ma(im, v.first, a.im);
              }
              if (beta != T(0)) {
                auto prior = V::load(y + out + r);
                re = V::ma(re, prior.first, beta);
                im = V::ma(im, prior.second, beta);
              }
              V::store(y + out + r, re, im);
            }
            for (; r < rhs; ++r) {
              Cx<T> sum{};
              for (uint32_t k = 0; k < 6; ++k) {
                auto a = blocks[(size_t(i) * 2 + ch) * 36 + q * 6 + k];
                sum = plus(
                    sum,
                    cm(Cx<T>{T(a.re), T(a.im)},
                       x[((size_t(i) * 3 + k % 3) * 4 + ch * 2 + k / 3) * rhs +
                         r]));
              }
              if (beta != T(0)) {
                sum.re += beta * y[out + r].re;
                sum.im += beta * y[out + r].im;
              }
              y[out + r] = sum;
            }
          }
    });
  }

  template <class T>
  void apply_cpu(const System &a, const Cx<T> *x, Cx<T> *y, int workers) {
    if constexpr (std::is_same_v<T, float>) {
      if (a.apply32) {
        a.apply32(x, y, workers);
        return;
      }
    } else {
      if (a.apply64) {
        a.apply64(x, y, workers);
        return;
      }
    }
    const auto &p = a.p;
    const Cx<T> *coeff;
    if constexpr (std::is_same_v<T, float>)
      coeff = a.af.data();
    else
      coeff = a.ad.data();
    using V = SIMD<T>;
    chunks(p.rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (uint32_t i = first; i < last; ++i)
        for (uint32_t c = 0; c < 3; ++c)
          for (uint32_t s = 0; s < 4; ++s) {
            uint32_t r = 0;
            for (; r + V::width <= p.rhs; r += V::width) {
              auto re = V::zero(), im = V::zero();
              for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e) {
                uint32_t si = p.perm[e * 4 + s];
                auto ph = p.phase[e * 4 + s];
                for (uint32_t k = 0; k < 3; ++k) {
                  auto av =
                      cm(coeff[e * 9 + c * 3 + k], Cx<T>{T(ph.re), T(ph.im)});
                  auto xv = V::load(
                      x + ((size_t(p.col[e]) * 3 + k) * 4 + si) * p.rhs + r);
                  re = V::ma(re, xv.first, av.re);
                  re = V::ms(re, xv.second, av.im);
                  im = V::ma(im, xv.second, av.re);
                  im = V::ma(im, xv.first, av.im);
                }
              }
              V::store(y + ((size_t(i) * 3 + c) * 4 + s) * p.rhs + r, re, im);
            }
            for (; r < p.rhs; ++r) {
              Cx<T> z{};
              for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e)
                for (uint32_t k = 0; k < 3; ++k) {
                  auto ph = p.phase[e * 4 + s];
                  auto av =
                      cm(coeff[e * 9 + c * 3 + k], Cx<T>{T(ph.re), T(ph.im)});
                  z = plus(z, cm(av, x[((size_t(p.col[e]) * 3 + k) * 4 +
                                        p.perm[e * 4 + s]) *
                                           p.rhs +
                                       r]));
                }
              y[((size_t(i) * 3 + c) * 4 + s) * p.rhs + r] = z;
            }
          }
    });
    if constexpr (std::is_same_v<T, float>) {
      if (!a.local32.empty())
        local_blocks(a.local32, x, y, p.rows, p.rhs, workers, T(1));
    } else {
      if (!a.local64.empty())
        local_blocks(a.local64, x, y, p.rows, p.rhs, workers, T(1));
    }
  }

  // Fixed unitary transports: A = mass*I + sum_mu(2I - T_mu - T_mu^*).
  // T_mu is a unitary periodic shift with dense color and spin-permutation
  // factors. Thus A is Hermitian positive definite and spectrum lies in
  // [mass,mass+16].
  System make_system(uint32_t l, uint32_t rhs, double mass, int workers,
                     bool zero_rhs = false) {
    System a;
    a.mass = mass;
    a.p = make_problem(l, rhs, true);
    auto &p = a.p;
    a.ad.resize(p.a.size());
    a.af.resize(p.a.size());
    using CD = std::complex<double>;
    const double pi = std::acos(-1.0);
    const CD roots[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (uint32_t i = 0; i < p.rows; ++i)
      for (uint32_t d = 0; d < 9; ++d) {
        uint32_t e = p.row[i] + d, mu = d ? (d - 1) / 2 : 0, power = mu + 1;
        bool forward = d % 2;
        for (uint32_t c = 0; c < 3; ++c)
          for (uint32_t k = 0; k < 3; ++k) {
            CD z = 0;
            if (!d)
              z = c == k ? 8 + mass : 0;
            else
              for (uint32_t j = 0; j < 3; ++j) {
                CD eig = j == 0   ? CD(1)
                         : j == 1 ? roots[power % 4]
                                  : roots[(2 * power) % 4];
                if (!forward)
                  eig = std::conj(eig);
                z -= eig * std::polar(1.0 / 3.0,
                                      2 * pi * (int(c) - int(k)) * j / 3.0);
              }
            size_t idx = e * 9 + c * 3 + k;
            a.ad[idx] = {z.real(), z.imag()};
            a.af[idx] = {float(z.real()), float(z.imag())};
            p.a[idx] = {float(z.real()), float(z.imag())};
          }
        for (uint32_t s = 0; s < 4; ++s) {
          p.perm[e * 4 + s] = !d ? s : (s + (forward ? power : 4 - power)) % 4;
          CD ph = !d ? CD(1) : (forward ? roots[mu] : std::conj(roots[mu]));
          p.phase[e * 4 + s] = {float(ph.real()), float(ph.imag())};
        }
      }
    // A known, independent solution with double-precision components.
    std::mt19937_64 rng(31771);
    std::uniform_real_distribution<double> dist(-1, 1);
    a.truth.resize(p.outputs());
    a.b.resize(p.outputs());
    for (auto &z : a.truth)
      z = {dist(rng), dist(rng)};
    if (zero_rhs && rhs > 1)
      for (size_t j = 0; j < p.rows * 12; ++j)
        a.truth[j * rhs + rhs - 1] = {0, 0}; // self-test zero RHS locking
    apply_cpu(a, a.truth.data(), a.b.data(), workers);
    a.bnorm.assign(rhs, 0);
    for (size_t j = 0; j < p.outputs(); ++j)
      a.bnorm[j % rhs] += a.b[j].re * a.b[j].re + a.b[j].im * a.b[j].im;
    for (auto &v : a.bnorm)
      v = std::sqrt(v);
    return a;
  }

  // Validate the mathematical assumptions independently of the CG recurrence.
  // Every directed transport must have the matching adjoint reverse edge, and
  // its color block must be unitary. These checks also work with duplicate
  // neighbors.
  void check_spd_structure(const System &a) {
    const auto &p = a.p;
    using CD = std::complex<double>;
    auto z = [](D v) { return CD(v.re, v.im); };
    for (uint32_t i = 0; i < p.rows; ++i)
      for (uint32_t d = 1; d < 9; ++d) {
        uint32_t e = p.row[i] + d,
                 back = p.row[p.col[e]] + (d % 2 ? d + 1 : d - 1);
        if (p.col[back] != i)
          throw std::runtime_error("reverse edge mismatch");
        for (uint32_t s = 0; s < 4; ++s) {
          uint32_t si = p.perm[e * 4 + s];
          if (p.perm[back * 4 + si] != s)
            throw std::runtime_error("spin inverse mismatch");
          auto ph = p.phase[e * 4 + s], pb = p.phase[back * 4 + si];
          for (uint32_t c = 0; c < 3; ++c)
            for (uint32_t k = 0; k < 3; ++k)
              if (std::abs(z(a.ad[e * 9 + c * 3 + k]) * CD(ph.re, ph.im) -
                           std::conj(z(a.ad[back * 9 + k * 3 + c]) *
                                     CD(pb.re, pb.im))) > 1e-12)
                throw std::runtime_error("operator is not Hermitian");
        }
        if (i == 0)
          for (uint32_t c = 0; c < 3; ++c)
            for (uint32_t k = 0; k < 3; ++k) {
              CD sum = 0;
              for (uint32_t j = 0; j < 3; ++j)
                sum += std::conj(z(a.ad[e * 9 + j * 3 + c])) *
                       z(a.ad[e * 9 + j * 3 + k]);
              if (std::abs(sum - CD(c == k ? 1 : 0)) > 1e-12)
                throw std::runtime_error("color transport is not unitary");
            }
      }
    CD energy = 0;
    double norm = 0;
    for (size_t j = 0; j < a.truth.size(); ++j) {
      energy += std::conj(z(a.truth[j])) * z(a.b[j]);
      norm += std::norm(z(a.truth[j]));
    }
    if (std::abs(energy.imag()) > 1e-12 * std::max(1.0, energy.real()) ||
        energy.real() < a.mass * norm * (1 - 1e-12))
      throw std::runtime_error("positive-definiteness energy check failed");
  }

  struct InnerInfo {
    int min_iterations = 0, max_iterations = 0, dispatched = 0;
    std::vector<int> rhs_iterations{};
  };
  template <class T> struct CPUCG {
    const System &a;
    int workers;
    std::vector<Cx<T>> delta, r, p, ap;
    std::vector<T> partial, rho, initial, pap, next, alpha, beta;
    std::vector<int> active, iterations;
    CPUCG(const System &a, int w)
        : a(a), workers(w), delta(a.p.outputs()), r(delta.size()),
          p(delta.size()), ap(delta.size()), partial(w * a.p.rhs), rho(a.p.rhs),
          initial(rho), pap(rho), next(rho), alpha(rho), beta(rho),
          active(a.p.rhs), iterations(a.p.rhs) {}
    void dot(const std::vector<Cx<T>> &v, const std::vector<Cx<T>> &w,
             std::vector<T> &out) {
      uint32_t rhs = a.p.rhs;
      chunks(a.p.rows * 12, workers,
             [&](int worker, uint32_t first, uint32_t last) {
               T *sums = partial.data() + worker * rhs;
               std::fill(sums, sums + rhs, T(0));
               for (size_t j = first; j < last; ++j)
                 for (uint32_t k = 0; k < rhs; ++k) {
                   auto z = v[j * rhs + k], q = w[j * rhs + k];
                   sums[k] += z.re * q.re + z.im * q.im;
                 }
             });
      for (uint32_t k = 0; k < rhs; ++k) {
        out[k] = 0;
        for (int w = 0; w < workers; ++w)
          out[k] += partial[w * rhs + k];
      }
    }
    InnerInfo run(const std::vector<D> &residual, double tol, int limit) {
      uint32_t rhs = a.p.rhs;
      chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
          r[j] = p[j] = {T(residual[j].re), T(residual[j].im)};
          delta[j] = {};
        }
      });
      dot(r, r, rho);
      initial = rho;
      std::fill(iterations.begin(), iterations.end(), 0);
      for (uint32_t k = 0; k < rhs; ++k)
        active[k] = rho[k] > 0;
      int it = 0;
      while (std::any_of(active.begin(), active.end(),
                         [](int v) { return v != 0; })) {
        if (it >= limit)
          throw std::runtime_error("CPU CG iteration limit");
        apply_cpu(a, p.data(), ap.data(), workers);
        dot(p, ap, pap);
        for (uint32_t k = 0; k < rhs; ++k) {
          if (active[k] && (!(pap[k] > 0) || !std::isfinite(pap[k])))
            throw std::runtime_error("CPU CG breakdown");
          alpha[k] = active[k] ? rho[k] / pap[k] : 0;
        }
        chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
            T t = alpha[j % rhs];
            delta[j].re += t * p[j].re;
            delta[j].im += t * p[j].im;
            r[j].re -= t * ap[j].re;
            r[j].im -= t * ap[j].im;
          }
        });
        dot(r, r, next);
        for (uint32_t k = 0; k < rhs; ++k) {
          if (active[k])
            ++iterations[k];
          if (!std::isfinite(next[k]))
            throw std::runtime_error("nonfinite CPU CG norm");
          if (next[k] <= T(tol * tol) * initial[k])
            active[k] = 0;
          beta[k] = active[k] ? next[k] / rho[k] : 0;
        }
        chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
            uint32_t k = j % rhs;
            p[j] = active[k] ? Cx<T>{r[j].re + beta[k] * p[j].re,
                                     r[j].im + beta[k] * p[j].im}
                             : Cx<T>{};
          }
        });
        rho.swap(next);
        ++it;
      }
      return {*std::min_element(iterations.begin(), iterations.end()),
              *std::max_element(iterations.begin(), iterations.end()), it};
    }
    const Cx<T> *correction() const { return delta.data(); }
  };

  struct GPUCG {
    GPUProblem matrix;
    MetalRun &m;
    std::array<id<MTLComputePipelineState>, 6> pipes;
    id<MTLBuffer> r, delta, partial, scalars;
    struct Params {
      uint32_t components, parts, limit;
      float tolerance2;
    } params;
    int chunk, workers;
    uint32_t rhs;
    GPUCG(MetalRun &m, const System &a, int chunk, int workers)
        : matrix(m, a.p), m(m), chunk(chunk), workers(workers), rhs(a.p.rhs) {
      params = {a.p.rows * 12, (a.p.rows * 12 + 255) / 256, 0, 0};
      auto values = [MTLFunctionConstantValues new];
      [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
      const char *names[] = {"cg_partial", "cg_init", "cg_alpha",
                             "cg_update",  "cg_beta", "cg_direction"};
      NSError *error = nil;
      for (int j = 0; j < 6; ++j) {
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
          throw std::runtime_error("GPU CG allocation failed");
        return b;
      };
      r = alloc(a.p.outputs() * 8);
      delta = alloc(a.p.outputs() * 8);
      partial = alloc(size_t(params.parts) * rhs * 4);
      scalars = alloc(rhs * 7 * 4);
    }
    InnerInfo run(const std::vector<D> &residual, double tol, int limit) {
      params.limit = limit;
      params.tolerance2 = float(tol * tol);
      auto *rr = static_cast<Cx<float> *>([r contents]);
      auto *pp = static_cast<Cx<float> *>([matrix.buffers[3] contents]);
      auto *dd = static_cast<Cx<float> *>([delta contents]);
      chunks(
          params.components, workers, [&](int, uint32_t first, uint32_t last) {
            for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
              rr[j] = pp[j] = {float(residual[j].re), float(residual[j].im)};
              dd[j] = {};
            }
          });
      auto *ss = static_cast<float *>([scalars contents]);
      int dispatched = 0;
      for (int offset = 0; offset < limit; offset += chunk) {
        @autoreleasepool {
          auto cmd = [m.queue commandBuffer];
          auto enc = [cmd computeCommandEncoder];
          auto dispatch = [&](int pipe, std::initializer_list<id<MTLBuffer>> bs,
                              size_t n) {
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
          if (offset == 0) {
            dispatch(0, {r, r, partial}, size_t(params.parts) * rhs);
            dispatch(1, {partial, scalars}, rhs);
          }
          for (int it = offset; it < std::min(limit, offset + chunk); ++it) {
            [enc setComputePipelineState:matrix.pipeline];
            for (size_t j = 0; j < matrix.buffers.size(); ++j)
              [enc setBuffer:matrix.buffers[j] offset:0 atIndex:j];
            [enc setBytes:&matrix.rows length:4 atIndex:7];
            [enc dispatchThreads:MTLSizeMake(matrix.outputs, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
            [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
            dispatch(0, {matrix.buffers[3], matrix.buffers[4], partial},
                     size_t(params.parts) * rhs);
            dispatch(2, {partial, scalars}, rhs);
            dispatch(3,
                     {delta, r, matrix.buffers[3], matrix.buffers[4], scalars},
                     matrix.outputs);
            dispatch(0, {r, r, partial}, size_t(params.parts) * rhs);
            dispatch(4, {partial, scalars}, rhs);
            dispatch(5, {r, matrix.buffers[3], scalars}, matrix.outputs);
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
          if (ss[6 * rhs + k])
            throw std::runtime_error("GPU CG breakdown");
          active |= ss[4 * rhs + k] > 0;
        }
        if (!active)
          break;
      }
      int lo = limit, hi = 0;
      for (uint32_t k = 0; k < rhs; ++k) {
        if (ss[k] > params.tolerance2 * ss[rhs + k])
          throw std::runtime_error("GPU CG iteration limit");
        lo = std::min(lo, int(ss[5 * rhs + k]));
        hi = std::max(hi, int(ss[5 * rhs + k]));
      }
      return {lo, hi, dispatched};
    }
    const Cx<float> *correction() const {
      return static_cast<Cx<float> *>([delta contents]);
    }
  };

  struct Trace {
    int outer;
    double max_residual;
    InnerInfo inner;
  };
  struct SolveResult {
    double total_ms = 0, inner_ms = 0, residual_ms = 0, update_ms = 0,
           max_residual = 0, solution_error = 0;
    int outer = 0, inner_max_sum = 0, dispatched = 0;
    std::vector<Trace> trace;
  };
  struct OuterWorkspace {
    std::vector<D> x, residual, ax;
    std::vector<double> sums, relative;
    OuterWorkspace(const System &a, int workers)
        : x(a.p.outputs()), residual(x.size()), ax(x.size()),
          sums(workers * a.p.rhs), relative(a.p.rhs) {}
  };

  template <class Inner>
  SolveResult solve(const System &a, Inner &inner, OuterWorkspace &w,
                    int workers, double tol, double inner_tol, int limit) {
    SolveResult out;
    auto start = Clock::now();
    uint32_t rhs = a.p.rhs;
    std::fill(w.x.begin(), w.x.end(), D{});
    w.residual = a.b;
    for (int outer = 0; outer < 12; ++outer) {
      auto t = Clock::now();
      InnerInfo info = inner.run(w.residual, inner_tol, limit);
      out.inner_ms += seconds(t) * 1000;
      t = Clock::now();
      auto *dx = inner.correction();
      chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j) {
          w.x[j].re += double(dx[j].re);
          w.x[j].im += double(dx[j].im);
        }
      });
      out.update_ms += seconds(t) * 1000;
      t = Clock::now();
      apply_cpu(a, w.x.data(), w.ax.data(), workers);
      chunks(a.p.rows * 12, workers,
             [&](int worker, uint32_t first, uint32_t last) {
               auto *sums = w.sums.data() + worker * rhs;
               std::fill(sums, sums + rhs, 0.0);
               for (size_t j = first; j < last; ++j)
                 for (uint32_t k = 0; k < rhs; ++k) {
                   size_t idx = j * rhs + k;
                   w.residual[idx] = {a.b[idx].re - w.ax[idx].re,
                                      a.b[idx].im - w.ax[idx].im};
                   sums[k] += w.residual[idx].re * w.residual[idx].re +
                              w.residual[idx].im * w.residual[idx].im;
                 }
             });
      out.max_residual = 0;
      for (uint32_t k = 0; k < rhs; ++k) {
        double sum = 0;
        for (int j = 0; j < workers; ++j)
          sum += w.sums[j * rhs + k];
        // For a zero RHS use absolute residual; the exact solution stays zero.
        double rel = std::sqrt(sum) / (a.bnorm[k] > 0 ? a.bnorm[k] : 1);
        if (!std::isfinite(rel))
          throw std::runtime_error("nonfinite true residual");
        w.relative[k] = rel;
        out.max_residual = std::max(out.max_residual, rel);
      }
      out.residual_ms += seconds(t) * 1000;
      out.outer = outer + 1;
      out.inner_max_sum += info.max_iterations;
      out.dispatched += info.dispatched;
      out.trace.push_back({outer + 1, out.max_residual, info});
      if (out.max_residual <= tol)
        break;
      // Lock RHS that already satisfy the true FP64 criterion.
      chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t j = first; j < last; ++j)
          for (uint32_t k = 0; k < rhs; ++k)
            if (w.relative[k] <= tol)
              w.residual[j * rhs + k] = {};
      });
      if (outer == 11)
        throw std::runtime_error(
            "iterative refinement failed to reach true tolerance");
    }
    out.total_ms = seconds(start) * 1000;
    // Independent known-solution error, outside solve timing. Residual above is
    // fully checked.
    long double err = 0, ref = 0;
    for (size_t j = 0; j < w.x.size(); ++j) {
      double re = w.x[j].re - a.truth[j].re, im = w.x[j].im - a.truth[j].im;
      err += re * re + im * im;
      ref += a.truth[j].re * a.truth[j].re + a.truth[j].im * a.truth[j].im;
    }
    out.solution_error = std::sqrt(double(err / std::max(ref, 1e-300L)));
    double error_limit = a.solution_error_limit > 0
                             ? a.solution_error_limit
                             : tol * (a.mass + 16) / a.mass * 2;
    if (!std::isfinite(out.solution_error) || out.solution_error > error_limit)
      throw std::runtime_error("known-solution check failed");
    return out;
  }

  std::vector<uint32_t> list_arg(const std::string &s, int max) {
    std::vector<uint32_t> v;
    std::stringstream in(s);
    std::string part;
    while (std::getline(in, part, ',')) {
      int n = std::stoi(part);
      if (n < 1 || n > max)
        throw std::runtime_error("list value out of range");
      v.push_back(n);
    }
    if (v.empty())
      throw std::runtime_error("empty list");
    return v;
  }

} // namespace MetalCloverEnv
