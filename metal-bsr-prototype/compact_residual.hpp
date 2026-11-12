#pragma once
#include "chroma_bridge.hpp"
#include "compact_clover.hpp"

namespace MetalCloverEnv {

  namespace compact_fp64 {
    using Vector = float64x2_t;

    struct ComplexVector {
      Vector re, im;
    };

    // Width 1 handles the final RHS without reading a second element.
    template <int Width> inline ComplexVector load(const D *x) {
      if constexpr (Width == 2) {
        auto z = SIMD<double>::load(x);
        return {z.first, z.second};
      } else {
        return {vsetq_lane_f64(x->re, vdupq_n_f64(0), 0),
                vsetq_lane_f64(x->im, vdupq_n_f64(0), 0)};
      }
    }

    template <int Width> inline void store(D *x, ComplexVector z) {
      if constexpr (Width == 2) {
        SIMD<double>::store(x, z.re, z.im);
      } else {
        *x = {vgetq_lane_f64(z.re, 0), vgetq_lane_f64(z.im, 0)};
      }
    }

    template <int Width>
    inline void site_residual(const CompactClover &c, const D *x, const D *b,
                              D *residual, uint32_t site, uint32_t r,
                              double *sums, uint32_t source_tag = 0,
                              bool zero_source = false) {
      const uint32_t rhs = c.w.a.p.rhs;
      std::array<ComplexVector, 12> output{};

      // Each (1 +/- gamma_mu) projector has two independent spin components.
      // Transport those two components and reconstruct the lower components.
      for (uint32_t direction = 0; direction < 8; ++direction) {
        uint32_t mu = direction / 2;
        uint32_t neighbor = c.neighbors[size_t(site) * 8 + direction];
        double sign = direction % 2 == 0 ? -1 : 1;
        std::array<ComplexVector, 6> projected;
        for (uint32_t color = 0; color < 3; ++color) {
          for (uint32_t spin = 0; spin < 2; ++spin) {
            uint32_t partner = c.perm[mu * 4 + spin];
            auto phase = c.phase[mu * 4 + spin];
            auto top = load<Width>(
                x + (size_t(neighbor) * 12 + color * 4 + spin) * rhs + r);
            auto lower = load<Width>(
                x + (size_t(neighbor) * 12 + color * 4 + partner) * rhs + r);
            auto &p = projected[color * 2 + spin];
            p.re = vfmsq_n_f64(vfmaq_n_f64(top.re, lower.re, sign * phase.re),
                               lower.im, sign * phase.im);
            p.im = vfmaq_n_f64(vfmaq_n_f64(top.im, lower.im, sign * phase.re),
                               lower.re, sign * phase.im);
          }
        }
        for (uint32_t color = 0; color < 3; ++color) {
          for (uint32_t spin = 0; spin < 2; ++spin) {
            ComplexVector transported{};
            for (uint32_t k = 0; k < 3; ++k) {
              auto gauge = c.gauge<double>(site, direction, color, k);
              const auto &p = projected[k * 2 + spin];
              transported.re = vfmsq_n_f64(
                  vfmaq_n_f64(transported.re, p.re, gauge.re), p.im, gauge.im);
              transported.im = vfmaq_n_f64(
                  vfmaq_n_f64(transported.im, p.im, gauge.re), p.re, gauge.im);
            }
            auto &top = output[color * 4 + spin];
            top.re = vaddq_f64(top.re, transported.re);
            top.im = vaddq_f64(top.im, transported.im);
            uint32_t partner = c.perm[mu * 4 + spin];
            auto phase = c.phase[mu * 4 + spin];
            auto &lower = output[color * 4 + partner];
            lower.re = vfmaq_n_f64(
                vfmaq_n_f64(lower.re, transported.re, sign * phase.re),
                transported.im, sign * phase.im);
            lower.im = vfmsq_n_f64(
                vfmaq_n_f64(lower.im, transported.im, sign * phase.re),
                transported.re, sign * phase.im);
          }
        }
      }

      for (uint32_t color = 0; color < 3; ++color) {
        for (uint32_t spin = 0; spin < 4; ++spin) {
          uint32_t chirality = spin / 2;
          uint32_t row = (spin % 2) * 3 + color;
          auto &out = output[color * 4 + spin];
          for (uint32_t k = 0; k < 6; ++k) {
            auto a = c.onsite_value(c.parity.parity[site], c.parity.index[site],
                                    chirality, row, k);
            auto v = load<Width>(
                x +
                (size_t(site) * 12 + (k % 3) * 4 + chirality * 2 + k / 3) *
                    rhs +
                r);
            out.re = vfmsq_n_f64(vfmaq_n_f64(out.re, v.re, a.re), v.im, a.im);
            out.im = vfmaq_n_f64(vfmaq_n_f64(out.im, v.im, a.re), v.re, a.im);
          }
        }
      }

      Vector norm = Width == 2 ? vld1q_f64(sums + r)
                               : vsetq_lane_f64(sums[r], vdupq_n_f64(0), 0);
      for (uint32_t component = 0; component < 12; ++component) {
        size_t q = (size_t(site) * 12 + component) * rhs + r;
        ComplexVector source{};
        if (b)
          source = load<Width>(b + q);
        else if (!zero_source) {
          auto z = chroma_bridge::source(q, source_tag);
          double scale = 1 / std::sqrt(double(c.w.a.p.outputs()));
          source.re = vsetq_lane_f64(z.re * scale, source.re, 0);
          source.im = vsetq_lane_f64(z.im * scale, source.im, 0);
        }
        ComplexVector z{vsubq_f64(source.re, output[component].re),
                        vsubq_f64(source.im, output[component].im)};
        store<Width>(residual + q, z);
        norm = vaddq_f64(
            norm, vaddq_f64(vmulq_f64(z.re, z.re), vmulq_f64(z.im, z.im)));
      }
      if constexpr (Width == 2)
        vst1q_f64(sums + r, norm);
      else
        sums[r] = vgetq_lane_f64(norm, 0);
    }

    // Evaluate the original full equation in FP64, including both clover
    // parities. Fuse b - M*x and its per-RHS norm into the stencil's final
    // stores.
    inline void residual(const CompactClover &c, const D *x, const D *b, D *r,
                         double *sums, int workers) {
      const uint32_t rhs = c.w.a.p.rhs;
      chunks(
          c.w.a.p.rows, workers,
          [&](int worker, uint32_t first, uint32_t last) {
            std::vector<double> local_sums(rhs, 0.);
            for (uint32_t site = first; site < last; ++site) {
              uint32_t column = 0;
              for (; column + 2 <= rhs; column += 2)
                site_residual<2>(c, x, b, r, site, column, local_sums.data());
              if (column < rhs)
                site_residual<1>(c, x, b, r, site, column, local_sums.data());
            }
            std::copy(local_sums.begin(), local_sums.end(),
                      sums + size_t(worker) * rhs);
          });
    }
  } // namespace compact_fp64

} // namespace MetalCloverEnv
