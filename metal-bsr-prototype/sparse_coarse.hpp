#pragma once
#include "fgmres.hpp"
#include <memory>

// Exact sparse Galerkin assembly through even sites. No global dense probing.
struct SparseCoarse {
  const TwoLevel &h;
  uint32_t bs, rows, n;
  int workers;
  double setup_ms = 0;
  std::vector<uint32_t> row, col, diagonal;
  std::vector<D> a64, di64;
  std::vector<Cx<float>> a32, di32;
  SparseCoarse(const TwoLevel &h, double mass, double gauge, double csw,
               bool apbc, int workers)
      : h(h), bs(2 * h.nvec), rows(h.aggregates), n(h.coarse),
        workers(workers) {
    auto start = Clock::now();
    CloverSystem c(h.extent, 1, mass, gauge, csw, apbc, workers);
    const auto &hop = c.parity.eo;
    std::vector<uint32_t> off{0}, ag;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> reverse(rows);
    std::vector<std::vector<uint32_t>> graph(rows);
    for (uint32_t a = 0; a < rows; ++a)
      graph[a].push_back(a);
    for (uint32_t e = 0; e < h.rows; ++e) {
      std::vector<uint32_t> v;
      for (uint32_t k = hop.p.row[e]; k < hop.p.row[e + 1]; ++k)
        v.push_back(h.aggregate[hop.p.col[k]]);
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      for (auto a : v) {
        reverse[a].push_back({e, uint32_t(ag.size())});
        ag.push_back(a);
        graph[a].insert(graph[a].end(), v.begin(), v.end());
      }
      off.push_back(uint32_t(ag.size()));
    }
    row.push_back(0);
    diagonal.resize(rows);
    for (uint32_t a = 0; a < rows; ++a) {
      auto &v = graph[a];
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      for (auto b : v) {
        if (a == b)
          diagonal[a] = uint32_t(col.size());
        col.push_back(b);
      }
      row.push_back(uint32_t(col.size()));
    }
    const size_t block = size_t(bs) * bs;
    a64.resize(col.size() * block);
    std::vector<D> t(size_t(ag.size()) * 12 * bs), u(t.size());
    chunks(h.rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (uint32_t e = first; e < last; ++e) {
        for (uint32_t k = hop.p.row[e]; k < hop.p.row[e + 1]; ++k) {
          uint32_t odd = hop.p.col[k], a = h.aggregate[odd];
          auto p = std::lower_bound(ag.begin() + off[e],
                                    ag.begin() + off[e + 1], a) -
                   ag.begin();
          auto *out = t.data() + size_t(p) * 12 * bs;
          for (uint32_t color = 0; color < 3; ++color)
            for (uint32_t s = 0; s < 4; ++s)
              for (uint32_t ic = 0; ic < 3; ++ic) {
                uint32_t si = hop.p.perm[k * 4 + s], q = color * 4 + s,
                         qi = ic * 4 + si;
                auto ph = hop.p.phase[k * 4 + s];
                auto z = cm(hop.ad[k * 9 + color * 3 + ic], D{ph.re, ph.im});
                for (uint32_t v = 0; v < h.nvec; ++v) {
                  auto term =
                      cm(z, h.basis[(size_t(odd) * 12 + qi) * h.nvec + v]);
                  auto &y = out[q * bs + (si / 2) * h.nvec + v];
                  y = plus(y, term);
                }
              }
        }
        for (uint32_t p = off[e]; p < off[e + 1]; ++p)
          for (uint32_t ch = 0; ch < 2; ++ch)
            for (uint32_t q = 0; q < 6; ++q)
              for (uint32_t v = 0; v < bs; ++v) {
                D sum{};
                for (uint32_t k = 0; k < 6; ++k)
                  sum = plus(
                      sum,
                      cm(c.inverseEven[(size_t(e) * 2 + ch) * 36 + q * 6 + k],
                         t[(size_t(p) * 12 + (k % 3) * 4 + ch * 2 + k / 3) *
                               bs +
                           v]));
                double sign = ch ? -1 : 1;
                u[(size_t(p) * 12 + (q % 3) * 4 + ch * 2 + q / 3) * bs + v] = {
                    sign * sum.re, sign * sum.im};
              }
      }
    });
    // Rows have disjoint ownership, so assembly needs no atomics or
    // volume-sized thread copies.
    chunks(rows, workers, [&](int, uint32_t first, uint32_t last) {
      std::vector<D> gram(block), ap(12 * bs);
      for (uint32_t a = first; a < last; ++a) {
        for (auto entry : reverse[a]) {
          auto e = entry.first, p = entry.second;
          for (uint32_t k = off[e]; k < off[e + 1]; ++k) {
            auto pos = std::lower_bound(col.begin() + row[a],
                                        col.begin() + row[a + 1], ag[k]) -
                       col.begin();
            auto *out = a64.data() + size_t(pos) * block;
            D one{1, 0}, zero{};
            cblas_zgemm(
                CblasRowMajor, CblasConjTrans, CblasNoTrans, int(bs), int(bs),
                12, reinterpret_cast<const __LAPACK_double_complex *>(&one),
                reinterpret_cast<const __LAPACK_double_complex *>(
                    t.data() + size_t(p) * 12 * bs),
                int(bs),
                reinterpret_cast<const __LAPACK_double_complex *>(
                    u.data() + size_t(k) * 12 * bs),
                int(bs),
                reinterpret_cast<const __LAPACK_double_complex *>(&zero),
                reinterpret_cast<__LAPACK_double_complex *>(gram.data()),
                int(bs));
            for (uint32_t i = 0; i < bs; ++i)
              for (uint32_t j = 0; j < bs; ++j) {
                double sign = i < h.nvec ? -1 : 1;
                auto z = gram[i * bs + j];
                out[i * bs + j].re += sign * z.re;
                out[i * bs + j].im += sign * z.im;
              }
          }
        }
        auto *out = a64.data() + size_t(diagonal[a]) * block;
        for (uint32_t p = h.offsets[a]; p < h.offsets[a + 1]; ++p) {
          uint32_t site = h.sites[p];
          std::fill(ap.begin(), ap.end(), D{});
          for (uint32_t color = 0; color < 3; ++color)
            for (uint32_t s = 0; s < 4; ++s) {
              uint32_t ch = s / 2, q = (s % 2) * 3 + color;
              for (uint32_t v = 0; v < h.nvec; ++v)
                for (uint32_t k = 0; k < 6; ++k) {
                  auto z = cm(c.oddA[(size_t(site) * 2 + ch) * 36 + q * 6 + k],
                              h.basis[(size_t(site) * 12 + (k % 3) * 4 +
                                       ch * 2 + k / 3) *
                                          h.nvec +
                                      v]);
                  ap[(color * 4 + s) * bs + ch * h.nvec + v] =
                      plus(ap[(color * 4 + s) * bs + ch * h.nvec + v], z);
                }
            }
          for (uint32_t q = 0; q < 12; ++q)
            for (uint32_t i = 0; i < h.nvec; ++i) {
              auto z = h.basis[(size_t(site) * 12 + q) * h.nvec + i];
              z.im = -z.im;
              uint32_t ci = (q % 4 / 2) * h.nvec + i;
              for (uint32_t j = 0; j < bs; ++j)
                out[ci * bs + j] =
                    plus(out[ci * bs + j], cm(z, ap[q * bs + j]));
            }
        }
      }
    });
    di64.resize(size_t(rows) * block);
    chunks(rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (uint32_t a = first; a < last; ++a) {
        std::vector<CD> d(block);
        for (size_t j = 0; j < block; ++j) {
          auto z = a64[size_t(diagonal[a]) * block + j];
          d[j] = {z.re, z.im};
        }
        auto inv = dense_inverse(d, bs);
        for (size_t j = 0; j < block; ++j)
          di64[size_t(a) * block + j] = {inv[j].real(), inv[j].imag()};
      }
    });
    a32 = block_cast<float>(a64);
    di32 = block_cast<float>(di64);
    setup_ms = seconds(start) * 1000;
    std::cout << "Sparse Galerkin: rows=" << rows << " block=" << bs
              << " nnz_blocks=" << col.size()
              << " FP32_MiB=" << a32.size() * 8. / 1048576
              << " assembly_ms=" << setup_ms << std::endl;
  }
  template <class T>
  void apply(const Cx<T> *x, Cx<T> *y, uint32_t rhs,
             bool diagonal_inverse = false) const {
    using V = SIMD<T>;
    const auto &coeff = [&]() -> const std::vector<Cx<T>> & {
      if constexpr (std::is_same_v<T, double>)
        return diagonal_inverse ? di64 : a64;
      else
        return diagonal_inverse ? di32 : a32;
    }();
    chunks(rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (uint32_t a = first; a < last; ++a)
        for (uint32_t i = 0; i < bs; ++i) {
          uint32_t r = 0;
          auto each = [&](const auto &f) {
            if (diagonal_inverse) {
              for (uint32_t j = 0; j < bs; ++j)
                f(coeff[(size_t(a) * bs + i) * bs + j], a * bs + j);
            } else
              for (uint32_t e = row[a]; e < row[a + 1]; ++e)
                for (uint32_t j = 0; j < bs; ++j)
                  f(coeff[(size_t(e) * bs + i) * bs + j], col[e] * bs + j);
          };
          for (; r + V::width <= rhs; r += V::width) {
            auto re = V::zero(), im = V::zero();
            each([&](Cx<T> z, uint32_t q) {
              auto v = V::load(x + size_t(q) * rhs + r);
              re = V::ma(re, v.first, z.re);
              re = V::ms(re, v.second, z.im);
              im = V::ma(im, v.second, z.re);
              im = V::ma(im, v.first, z.im);
            });
            V::store(y + size_t(a * bs + i) * rhs + r, re, im);
          }
          for (; r < rhs; ++r) {
            Cx<T> sum{};
            each([&](Cx<T> z, uint32_t q) {
              sum = plus(sum, cm(z, x[size_t(q) * rhs + r]));
            });
            y[size_t(a * bs + i) * rhs + r] = sum;
          }
        }
    });
  }
};

template <class Matrix> struct CPUBlockCoarseSolve {
  const Matrix &a;
  System fake;
  std::unique_ptr<CPUBiCGStab<float>> bi;
  uint32_t rhs;
  int limit;
  double tol;
  std::vector<Cx<float>> b, scaled, tmp;
  uint64_t calls = 0, applications = 0, iterations = 0;
  double worst_true = 0;
  CPUBlockCoarseSolve(const Matrix &a, uint32_t rhs, int workers, int limit,
                      double tol)
      : a(a), rhs(rhs), limit(limit), tol(tol), b(size_t(a.n) * rhs),
        scaled(b.size()), tmp(b.size()) {
    fake.p.rows = a.n;
    fake.p.rhs = rhs;
    fake.p.colors = 1;
    fake.p.spins = 1;
    bi = std::make_unique<CPUBiCGStab<float>>(fake, workers);
    bi->allow_inexact = true;
    bi->apply_override = [this](const Cx<float> *x, Cx<float> *y) {
      this->a.apply(x, tmp.data(), this->rhs);
      this->a.apply(tmp.data(), y, this->rhs, true);
      ++applications;
    };
  }
  void reset() {
    calls = applications = iterations = 0;
    worst_true = 0;
  }
  void solve(const Cx<float> *in, Cx<float> *out) {
    ++calls;
    std::copy(in, in + b.size(), b.begin());
    a.apply(b.data(), scaled.data(), rhs, true);
    auto info = bi->run(scaled, tol, limit);
    iterations += info.max_iterations;
    std::copy(bi->delta.begin(), bi->delta.end(), out);
    a.apply(out, tmp.data(), rhs);
    ++applications;
    std::vector<double> e(rhs), n(rhs);
    for (size_t i = 0; i < b.size(); ++i) {
      auto z = minus_z(b[i], tmp[i]);
      e[i % rhs] += double(z.re) * z.re + double(z.im) * z.im;
      n[i % rhs] += double(b[i].re) * b[i].re + double(b[i].im) * b[i].im;
    }
    for (uint32_t r = 0; r < rhs; ++r) {
      double v = std::sqrt(e[r] / std::max(n[r], 1e-300));
      if (!std::isfinite(v))
        throw std::runtime_error("nonfinite coarse true residual");
      worst_true = std::max(worst_true, v);
    }
  }
};
using CPUCoarseSolve = CPUBlockCoarseSolve<SparseCoarse>;

template <class Matrix> struct GPUBlockCoarseSolve {
  const Matrix &a;
  MetalRun &m;
  uint32_t rhs;
  int limit;
  double tol;
  id<MTLBuffer> rows, cols, coeff, diag, r, shadow, p, v, t, delta, tmp, scaled,
      partial, second, q;
  std::array<id<MTLComputePipelineState>, 10> pipes;
  struct BP {
    uint32_t components, parts, limit;
    float tolerance2;
    uint32_t stage;
  } bp;
  struct SP {
    uint32_t n, bs, parts, mode;
  } sp;
  uint64_t calls = 0, applications = 0, iterations = 0, padded_iterations = 0;
  double worst_true = 0;
  GPUBlockCoarseSolve(const Matrix &a, MetalRun &m, uint32_t rhs, int limit,
                      double tol)
      : a(a), m(m), rhs(rhs), limit(limit), tol(tol) {
    auto put = [&](const void *data, size_t bytes) {
      auto b = [m.device newBufferWithLength:std::max(size_t(8), bytes)
                                     options:MTLResourceStorageModeShared];
      if (!b)
        throw std::runtime_error("sparse coarse Metal allocation failed");
      if (data)
        std::memcpy([b contents], data, bytes);
      return b;
    };
    rows = put(a.row.data(), a.row.size() * 4);
    cols = put(a.col.data(), a.col.size() * 4);
    coeff = put(a.a32.data(), a.a32.size() * 8);
    diag = put(a.di32.data(), a.di32.size() * 8);
    size_t bytes = size_t(a.n) * rhs * 8;
    r = put(nullptr, bytes);
    shadow = put(nullptr, bytes);
    p = put(nullptr, bytes);
    v = put(nullptr, bytes);
    t = put(nullptr, bytes);
    delta = put(nullptr, bytes);
    tmp = put(nullptr, bytes);
    scaled = put(nullptr, bytes);
    bp = {a.n, (a.n + 255) / 256, uint32_t(limit), float(tol * tol), 0};
    sp = {a.n, a.bs, bp.parts, 0};
    partial = put(nullptr, size_t(bp.parts) * rhs * 8);
    second = put(nullptr, size_t(bp.parts) * rhs * 8);
    q = put(nullptr, size_t(rhs) * 11 * 8);
    const char *names[] = {"bi_initial",   "bi_partial",    "bi_scalar",
                           "bi_direction", "bi_update",     "sc_apply",
                           "sc_diagonal",  "sc_difference", "sc_diagnostics",
                           "sc_copy"};
    auto values = [MTLFunctionConstantValues new];
    [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
    for (int j = 0; j < 10; ++j) {
      NSError *error = nil;
      auto f = [m.library
          newFunctionWithName:[NSString stringWithUTF8String:names[j]]
               constantValues:values
                        error:&error];
      if (!f)
        throw std::runtime_error([[error localizedDescription] UTF8String]);
      pipes[j] = [m.device newComputePipelineStateWithFunction:f error:&error];
      if (!pipes[j])
        throw std::runtime_error([[error localizedDescription] UTF8String]);
    }
  }
  void reset() {
    calls = applications = iterations = padded_iterations = 0;
    worst_true = 0;
  }
  void dispatch(id<MTLComputeCommandEncoder> enc, int pipe,
                std::initializer_list<id<MTLBuffer>> bs, size_t n,
                uint32_t mode = 0) {
    [enc setComputePipelineState:pipes[pipe]];
    size_t j = 0;
    for (auto b : bs) {
      [enc setBuffer:b offset:0 atIndex:j];
      ++j;
    }
    bp.stage = sp.mode = mode;
    if (pipe < 5)
      [enc setBytes:&bp length:sizeof(bp) atIndex:j];
    else
      [enc setBytes:&sp length:sizeof(sp) atIndex:j];
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
  }
  void apply(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x, id<MTLBuffer> y,
             bool preconditioned = true) {
    dispatch(enc, 5, {rows, cols, coeff, x, preconditioned ? tmp : y, q},
             size_t(a.n) * rhs, preconditioned ? 1 : 0);
    ++applications;
    if (preconditioned)
      dispatch(enc, 6, {diag, tmp, y}, size_t(a.n) * rhs);
  }
  void encode_product(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
                      id<MTLBuffer> y) {
    apply(enc, x, y, false);
  }
  void difference(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> b,
                  id<MTLBuffer> ax, id<MTLBuffer> out) {
    dispatch(enc, 7, {b, ax}, size_t(a.n) * rhs);
    if (out != ax)
      dispatch(enc, 9, {ax, out}, size_t(a.n) * rhs);
  }
  void encode(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> input,
              id<MTLBuffer> out) {
    ++calls;
    size_t n = size_t(a.n) * rhs, np = size_t(bp.parts) * rhs;
    auto dot = [&](id<MTLBuffer> x, id<MTLBuffer> y, id<MTLBuffer> result) {
      dispatch(enc, 1, {x, y, result}, np);
    };
    auto scalar = [&](uint32_t stage) {
      dispatch(enc, 2, {partial, second, q}, rhs, stage);
    };
    dispatch(enc, 6, {diag, input, scaled}, n);
    dispatch(enc, 0, {scaled, r, shadow, p, v, delta}, n);
    dot(r, r, partial);
    scalar(0);
    for (int i = 0; i < limit; ++i) {
      dot(shadow, r, partial);
      scalar(1);
      dispatch(enc, 3, {r, p, v, q}, n);
      apply(enc, p, v);
      dot(shadow, v, partial);
      scalar(2);
      dispatch(enc, 4, {delta, r, p, v, t, q}, n, 0);
      dot(r, r, partial);
      scalar(3);
      apply(enc, r, t);
      dot(t, r, partial);
      dot(t, t, second);
      scalar(4);
      dispatch(enc, 4, {delta, r, p, v, t, q}, n, 1);
      dot(r, r, partial);
      scalar(5);
      ++padded_iterations;
    }
    apply(enc, delta, tmp, false);
    dispatch(enc, 7, {input, tmp}, n);
    dot(tmp, tmp, partial);
    dot(input, input, second);
    dispatch(enc, 8, {partial, second, q}, rhs);
    dispatch(enc, 9, {delta, out}, n);
  }
  void check() {
    auto *s = static_cast<Cx<float> *>([q contents]);
    int hi = 0;
    for (uint32_t r = 0; r < rhs; ++r) {
      if (s[8 * rhs + r].re)
        throw std::runtime_error("coarse GPU BiCGStab breakdown");
      hi = std::max(hi, int(s[7 * rhs + r].re));
      auto z = s[10 * rhs + r];
      double ratio = std::sqrt(double(z.re) / std::max(double(z.im), 1e-300));
      if (!std::isfinite(ratio))
        throw std::runtime_error("nonfinite GPU coarse true residual");
      worst_true = std::max(worst_true, ratio);
    }
    iterations += hi;
  }
};
using GPUCoarseSolve = GPUBlockCoarseSolve<SparseCoarse>;
