#pragma once
#include "compact_solve.hpp"
#include "mapped_field.hpp"
#include "production_bridge.hpp"
#include <atomic>

namespace MetalCloverEnv {

  inline void pack_hermitian6(const M6 &a, double *out) {
    for (int i = 0; i < 6; ++i)
      out[i] = a[i * 6 + i].real();
    for (int i = 1; i < 6; ++i)
      for (int j = 0; j < i; ++j) {
        int p = 6 + i * (i - 1) + 2 * j;
        out[p] = a[i * 6 + j].real();
        out[p + 1] = a[i * 6 + j].imag();
      }
  }
  inline void pack_hermitian6(const M6 &a, float *out) {
    double tmp[36];
    pack_hermitian6(a, tmp);
    for (int i = 0; i < 36; ++i)
      out[i] = float(tmp[i]);
  }

  // Hermitian local blocks are HPD on this action, as checked during setup.
  // Unit-lower L and real D avoid square roots; U=L-dagger is implicit.
  inline M6 factor_ldl6(const M6 &a, double &minimum) {
    M6 f{};
    for (int i = 0; i < 6; ++i) {
      double d = a[i * 6 + i].real();
      for (int k = 0; k < i; ++k)
        d -= std::norm(f[i * 6 + k]) * f[k * 6 + k].real();
      if (!(d > 0) || !std::isfinite(d))
        throw std::runtime_error("LDL local pivot is not positive");
      minimum = std::min(minimum, d);
      f[i * 6 + i] = d;
      for (int j = i + 1; j < 6; ++j) {
        CD z = a[j * 6 + i];
        for (int k = 0; k < i; ++k)
          z -= f[j * 6 + k] * f[k * 6 + k].real() * std::conj(f[i * 6 + k]);
        f[j * 6 + i] = z / d;
      }
    }
    return f;
  }

  struct LowMemoryClover : CompactClover {

    MappedField &bridge;
    production::Header header;
    size_t gauge_offset;
    void raw_onsite(uint32_t site, D *out) const {
      std::array<CD, 144> cl{};
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = mu + 1; nu < 4; ++nu) {
          int m = mu + 1, n = nu + 1;
          auto aa = gauge_loop(w, site, {m, n, -m, -n}),
               b = gauge_loop(w, site, {n, -m, -n, m}),
               cc = gauge_loop(w, site, {-m, -n, m, n}),
               d = gauge_loop(w, site, {-n, m, n, -m});
          M3 q;
          for (int j = 0; j < 9; ++j)
            q[j] = aa[j] + b[j] + cc[j] + d[j];
          auto adj = m3_adj(q);
          for (int j = 0; j < 9; ++j)
            q[j] = (q[j] - adj[j]) * .125;
          for (int ss = 0; ss < 4; ++ss)
            for (int t = 0; t < 4; ++t) {
              CD spin = 0;
              for (int k = 0; k < 4; ++k)
                spin += w.gamma[mu][ss * 4 + k] * w.gamma[nu][k * 4 + t];
              for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                  cl[(ss * 3 + i) * 12 + t * 3 + j] +=
                      -.5 * (nu == 3 ? electric : magnetic) * spin *
                      q[i * 3 + j];
            }
        }
      for (int ch = 0; ch < 2; ++ch)
        for (int i = 0; i < 6; ++i)
          for (int j = 0; j < 6; ++j) {
            auto z = cl[(ch * 6 + i) * 12 + ch * 6 + j] +
                     CD(i == j ? header.mass + hopping[0] + hopping[1] +
                                     hopping[2] + 1
                               : 0);
            out[ch * 36 + i * 6 + j] = {z.real(), z.imag()};
          }
    }
    LowMemoryClover(MappedField &bridge, const production::Header &h,
                    size_t gauge_offset, int threads, bool antiperiodic = true)
        : bridge(bridge), header(h), gauge_offset(gauge_offset) {
      if (gauge_offset % alignof(M3))
        throw std::runtime_error("Use an aligned CHRMG003 bridge");
      auto start = Clock::now();
      workers = threads;
      packed = true;
      factorized = true;
      w.extent = h.shape[0];
      w.dimensions = h.shape;
      w.antiperiodic = antiperiodic;
      w.strength = .25;
      w.mapped_links = reinterpret_cast<const M3 *>(
          static_cast<const char *>(bridge.data) + gauge_offset);
      hopping = {h.nu / h.xi0, h.nu / h.xi0, h.nu / h.xi0, 1};
      magnetic = h.clover_r / h.xi0;
      electric = h.clover_t;
      for (auto &g : w.gamma)
        for (auto &z : g)
          z = std::conj(z);
      auto &a = w.a;
      auto &p = a.p;
      p.rows = h.volume();
      p.colors = 3;
      p.spins = 4;
      p.rhs = 1;
      a.mass = h.mass;
      for (uint32_t mu = 0; mu < 4; ++mu)
        for (uint32_t s = 0; s < 4; ++s)
          for (uint32_t k = 0; k < 4; ++k)
            if (std::abs(w.gamma[mu][s * 4 + k]) > 0) {
              perm[mu * 4 + s] = k;
              auto z = w.gamma[mu][s * 4 + k];
              phase[mu * 4 + s] = {z.real(), z.imag()};
            }
      parity.index.resize(p.rows);
      parity.parity.resize(p.rows);
      neighbors.resize(size_t(p.rows) * 8);
      for (uint32_t i = 0; i < p.rows; ++i) {
        uint32_t q = i, pa = 0;
        for (int mu = 0; mu < 4; ++mu) {
          pa += q % h.shape[mu];
          q /= h.shape[mu];
        }
        pa %= 2;
        parity.parity[i] = pa;
        auto &list = pa ? parity.odd : parity.even;
        parity.index[i] = uint32_t(list.size());
        list.push_back(i);
        for (uint32_t d = 0; d < 8; ++d)
          neighbors[size_t(i) * 8 + d] = w.neighbor(i, d / 2, d % 2 == 0);
      }
      uint32_t half = p.rows / 2;
      links32.resize(size_t(p.rows) * 36);
      oddA32.resize(size_t(half) * 36);
      inverseEven32.resize(oddA32.size());

      static std::atomic<uint64_t> serial{0};
      generation = ++serial;
      generated_clover = [this](uint32_t site, D *out) {
        raw_onsite(site, out);
      };
      std::vector<double> pivots(workers, 1e100);
      std::vector<std::string> errors(workers);
      uint32_t slab = w.stride(3) * 8;
      for (uint32_t begin = 0; begin < p.rows; begin += slab) {
        uint32_t count = std::min(slab, p.rows - begin);
        chunks(count, workers, [&](int worker, uint32_t first, uint32_t last) {
          try {
            for (uint32_t site = begin + first; site < begin + last; ++site) {
              for (int mu = 0; mu < 4; ++mu)
                for (int j = 0; j < 9; ++j) {
                  auto z = w.link_at(size_t(site) * 4 + mu)[j];
                  links32[(size_t(site) * 4 + mu) * 9 + j] = {
                      float(hopping[mu] * z.real()),
                      float(hopping[mu] * z.imag())};
                }
              std::array<D, 72> raw;
              raw_onsite(site, raw.data());
              bool odd = parity.parity[site];
              uint32_t ix = parity.index[site];
              for (int ch = 0; ch < 2; ++ch) {
                M6 block;
                for (int k = 0; k < 36; ++k) {
                  auto z = raw[ch * 36 + k];
                  block[k] = {z.re, z.im};
                }
                auto factors = factor_ldl6(block, pivots[worker]);
                auto &storage = odd ? oddA32 : inverseEven32;
                pack_hermitian6(factors,
                                reinterpret_cast<float *>(storage.data()) +
                                    (size_t(ix) * 2 + ch) * 36);
              }
            }
          } catch (const std::exception &e) {
            errors[worker] = e.what();
          }
        });
        for (const auto &error : errors)
          if (!error.empty())
            throw std::runtime_error(error);
        bridge.drop();
      }
      min_pivot = *std::min_element(pivots.begin(), pivots.end());
      for (auto *halfsystem : {&parity.eo, &parity.oe}) {
        halfsystem->mass = h.mass;
        halfsystem->p.rows = half;
        halfsystem->p.colors = 3;
        halfsystem->p.spins = 4;
        halfsystem->p.rhs = 1;
      }
      a.apply64 = [this](const D *x, D *y, int n) { apply(x, y, n, -1, true); };
      parity.eo.apply64 = [this](const D *x, D *y, int n) {
        apply(x, y, n, 0, false);
      };
      parity.oe.apply64 = [this](const D *x, D *y, int n) {
        apply(x, y, n, 1, false);
      };
      setup_ms = seconds(start) * 1000;
      std::cout << "LDL_CLOVER seconds=" << setup_ms / 1000
                << " CPU_FP64_clover=recomputed_from_original_FP64_gauge min_D="
                << min_pivot << std::endl;
    }
    void drop_coefficients() { bridge.drop(); }
  };

  inline double streamed_residual(LowMemoryClover &c, MappedField &x,
                                  MappedField &residual, uint32_t tag,
                                  int workers, bool probe = false,
                                  Cx<float> *next_residual = nullptr) {
    const D *values = static_cast<const D *>(x.data);
    D *out = static_cast<D *>(residual.data);
    uint32_t slab = c.w.stride(3) * 8;
    double sum = 0, reference = 0;
    const D *expected = reinterpret_cast<const D *>(
        static_cast<const char *>(c.bridge.data) + c.gauge_offset +
        size_t(c.w.a.p.rows) * 36 * 16);
    for (uint32_t begin = 0; begin < c.w.a.p.rows; begin += slab) {
      uint32_t count = std::min(slab, c.w.a.p.rows - begin);
      std::vector<double> sums(workers), refs(workers);
      chunks(count, workers, [&](int worker, uint32_t first, uint32_t last) {
        double n = 0, ref = 0, error = 0;
        for (uint32_t site = begin + first; site < begin + last; ++site) {
          compact_fp64::site_residual<1>(c, values, nullptr, out, site, 0, &n,
                                         tag, probe);
          if (probe)
            for (int q = 0; q < 12; ++q) {
              size_t i = size_t(site) * 12 + q;
              auto z = expected[i], r = out[i];
              ref += z.re * z.re + z.im * z.im;
              error +=
                  (r.re + z.re) * (r.re + z.re) + (r.im + z.im) * (r.im + z.im);
            }
        }
        sums[worker] = probe ? error : n;
        refs[worker] = ref;
      });
      for (auto v : sums)
        sum += v;
      for (auto v : refs)
        reference += v;
      if (next_residual)
        chunks(count * 12, workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = size_t(begin) * 12 + first;
               i < size_t(begin) * 12 + last; ++i)
            next_residual[i] = {float(out[i].re), float(out[i].im)};
        });
      residual.drop(size_t(begin) * 12 * 16, size_t(count) * 12 * 16, true);
      x.drop();
      c.drop_coefficients();
    }
    if (probe)
      return std::sqrt(std::max(0., sum) / reference);
    return std::sqrt(sum);
  }

} // namespace MetalCloverEnv
