#pragma once
#include "chroma_bridge.hpp"
#include "mg_real_tune.hpp"

namespace MetalCloverEnv {

  // Like mgproton's homogeneous solve construction: z ~= S^{-1}(S y),
  // v = y-z. The benchmark RHS never participates in basis setup.
  struct NearNullOptions {
    int candidates = 24, iterations = 128, restart = 12, warm = 32, batch = 0;
    double tolerance = .001;
  };

  struct BasisSVDStats {
    double retained_energy = 0, minimum_retained_energy = 1;
    double input_energy = 0, kept_energy = 0;
    double minimum_rank_ratio = 1, maximum_orthogonality_error = 0;
  };

  // Candidate storage is FP32, component-major. Each aggregate/chirality is
  // compressed independently by a CPU FP64 thin SVD, retaining leading LEFT
  // singular vectors. Full candidate normalization precedes this local SVD.
  inline BasisSVDStats
  compress_near_null_svd(Transfer &h, std::vector<Cx<float>> &candidates,
                         uint32_t count, std::ostream *spectrum = nullptr,
                         uint32_t stride = 0, uint32_t carry = 0) {
    if (!stride)
      stride = count;
    uint32_t keep = carry ? carry : h.nvec;
    if (count < keep || stride < count ||
        candidates.size() != size_t(h.rows) * h.inbs * stride)
      throw std::runtime_error("invalid near-null candidate dimensions");
    if (!carry)
      h.basis.assign(size_t(h.rows) * h.inbs * h.nvec, {});
    std::vector<double> singular(size_t(h.aggregates) * 2 * count);
    std::vector<double> total(h.workers), kept(h.workers),
        minimum(h.workers, 1), rank(h.workers, 1), orthogonality(h.workers);
    std::vector<std::string> errors(h.workers);
    chunks(h.aggregates, h.workers,
           [&](int worker, uint32_t first, uint32_t last) {
             try {
               for (uint32_t ag = first; ag < last; ++ag)
                 for (uint32_t ch = 0; ch < 2; ++ch) {
                   std::vector<size_t> components;
                   for (uint32_t p = h.offsets[ag]; p < h.offsets[ag + 1]; ++p)
                     for (uint32_t q = 0; q < h.inbs; ++q)
                       if (h.chirality(q) == ch)
                         components.push_back(size_t(h.sites[p]) * h.inbs + q);
                   __LAPACK_int m = __LAPACK_int(components.size()), n = count,
                                lda = m, ldu = m, ldvt = 1, info = 0,
                                lwork = -1;
                   if (m < n)
                     throw std::runtime_error(
                         "more SVD candidates than local rows");
                   std::vector<__LAPACK_double_complex> a(size_t(m) * n),
                       u(a.size());
                   std::vector<double> s(n), rwork(5 * n);
                   for (int j = 0; j < n; ++j)
                     for (int i = 0; i < m; ++i) {
                       auto z = candidates[components[i] * stride + j];
                       a[i + size_t(m) * j] = CD(z.re, z.im);
                     }
                   char jobu = 'S', jobvt = 'N';
                   __LAPACK_double_complex query, unused;
                   zgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s.data(),
                           u.data(), &ldu, &unused, &ldvt, &query, &lwork,
                           rwork.data(), &info);
                   if (info)
                     throw std::runtime_error("SVD workspace query failed");
                   lwork = __LAPACK_int(std::ceil(CD(query).real()));
                   std::vector<__LAPACK_double_complex> work(lwork);
                   zgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s.data(),
                           u.data(), &ldu, &unused, &ldvt, work.data(), &lwork,
                           rwork.data(), &info);
                   if (info || !std::isfinite(s[0]) || !(s[0] > 0) ||
                       !(s[keep - 1] > 1e-12 * s[0]))
                     throw std::runtime_error("rank-deficient near-null SVD");
                   double sum = 0, retained = 0;
                   for (int j = 0; j < n; ++j) {
                     singular[(size_t(ag) * 2 + ch) * count + j] = s[j];
                     sum += s[j] * s[j];
                     if (j < int(keep))
                       retained += s[j] * s[j];
                   }
                   total[worker] += sum;
                   kept[worker] += retained;
                   minimum[worker] = std::min(minimum[worker], retained / sum);
                   rank[worker] = std::min(rank[worker], s[keep - 1] / s[0]);
                   for (uint32_t j = 0; j < keep; ++j)
                     for (int i = 0; i < m; ++i) {
                       CD z = u[i + size_t(m) * j];
                       if (carry)
                         candidates[components[i] * stride + j] = {
                             float(z.real() * s[j]), float(z.imag() * s[j])};
                       else
                         h.basis[components[i] * h.nvec + j] = {z.real(),
                                                                z.imag()};
                     }
                   for (uint32_t j = 0; j < keep; ++j)
                     for (uint32_t k = 0; k <= j; ++k) {
                       CD dot = 0;
                       for (int i = 0; i < m; ++i)
                         dot += std::conj(CD(u[i + size_t(m) * k])) *
                                CD(u[i + size_t(m) * j]);
                       orthogonality[worker] =
                           std::max(orthogonality[worker],
                                    std::abs(dot - CD(j == k ? 1 : 0)));
                     }
                 }
             } catch (const std::exception &e) {
               errors[worker] = e.what();
             }
           });
    for (const auto &e : errors)
      if (!e.empty())
        throw std::runtime_error(e);
    BasisSVDStats stats;
    double sum = 0, retained = 0;
    for (int w = 0; w < h.workers; ++w) {
      sum += total[w];
      retained += kept[w];
      stats.minimum_retained_energy =
          std::min(stats.minimum_retained_energy, minimum[w]);
      stats.minimum_rank_ratio = std::min(stats.minimum_rank_ratio, rank[w]);
      stats.maximum_orthogonality_error =
          std::max(stats.maximum_orthogonality_error, orthogonality[w]);
    }
    stats.input_energy = sum;
    stats.kept_energy = retained;
    stats.retained_energy = retained / sum;
    if (stats.maximum_orthogonality_error > 1e-12)
      throw std::runtime_error("near-null prolongator is not orthogonal");
    if (!carry)
      h.basis32 = block_cast<float>(h.basis);
    if (spectrum) {
      *spectrum << "aggregate,chirality,index,singular_value\n";
      for (uint32_t ag = 0; ag < h.aggregates; ++ag)
        for (uint32_t ch = 0; ch < 2; ++ch)
          for (uint32_t j = 0; j < count; ++j)
            *spectrum << ag << ',' << ch << ',' << j << ','
                      << singular[(size_t(ag) * 2 + ch) * count + j] << '\n';
    }
    return stats;
  }

  inline void check_near_null_svd() {
    Transfer h(4, 4, 2, LatticeShape{2, 2, 1, 1}, LatticeShape{2, 2, 1, 1},
               false, 2);
    std::vector<Cx<float>> candidates(4 * 4 * 3);
    for (uint32_t ch = 0; ch < 2; ++ch) {
      candidates[(ch * 2) * 3] = {0, 3};
      candidates[(ch * 2 + 1) * 3 + 1] = {1.2f, 1.6f};
      candidates[(4 + ch * 2) * 3 + 2] = {1, 0};
    }
    auto s = compress_near_null_svd(h, candidates, 3);
    if (std::abs(s.retained_energy - 13. / 14.) > 1e-8)
      throw std::runtime_error("analytic SVD retained-energy test failed");
    for (uint32_t ch = 0; ch < 2; ++ch)
      for (uint32_t j = 0; j < 3; ++j) {
        double residual = 0;
        CD coefficients[2]{};
        for (uint32_t v = 0; v < 2; ++v)
          for (uint32_t q = 0; q < 16; ++q)
            if (h.chirality(q % 4) == ch) {
              auto a = h.basis[q * 2 + v];
              auto b = candidates[q * 3 + j];
              coefficients[v] += std::conj(CD(a.re, a.im)) * CD(b.re, b.im);
            }
        for (uint32_t q = 0; q < 16; ++q)
          if (h.chirality(q % 4) == ch) {
            auto b = candidates[q * 3 + j];
            CD error(b.re, b.im);
            for (uint32_t v = 0; v < 2; ++v) {
              auto a = h.basis[q * 2 + v];
              error -= CD(a.re, a.im) * coefficients[v];
            }
            residual += std::norm(error);
          }
        if (std::abs(residual - (j == 2 ? 1. : 0.)) > 1e-12)
          throw std::runtime_error("analytic SVD projector test failed");
      }
    std::fill(candidates.begin(), candidates.end(), Cx<float>{});
    bool rejected = false;
    try {
      compress_near_null_svd(h, candidates, 3);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    if (!rejected)
      throw std::runtime_error("rank-deficient SVD was accepted");
    {
      Transfer test(4, 4, 2, LatticeShape{2, 2, 1, 1}, LatticeShape{2, 2, 1, 1},
                    false, 2);
      std::vector<Cx<float>> pool(16 * 4);
      for (uint32_t ch = 0; ch < 2; ++ch)
        for (uint32_t j = 0; j < 4; ++j)
          pool[((j / 2) * 4 + ch * 2 + j % 2) * 4 + j] = {float(4 - j), 0};
      auto carry = compress_near_null_svd(test, pool, 4, nullptr, 4, 3);
      if (!test.basis.empty() || std::abs(carry.kept_energy - 58) > 1e-10)
        throw std::runtime_error("weighted SVD carry test failed");
      for (uint32_t ch = 0; ch < 2; ++ch)
        pool[(4 + ch * 2 + 1) * 4 + 3] = {0, 5};
      auto final = compress_near_null_svd(test, pool, 4);
      if (std::abs(final.kept_energy - 82) > 1e-5 ||
          std::abs(final.input_energy - 108) > 1e-5)
        throw std::runtime_error("incremental weighted SVD energy test failed");
    }
    std::cout << "PASS analytic complex/chiral SVD truncation, projector, "
                 "orthogonality and rank rejection"
              << std::endl;
  }

  template <class Schur>
  void build_gpu_near_null_svd(Transfer &h, Schur &op, const NearNullOptions &o,
                               std::ostream &vectors_log,
                               std::ostream &spectrum) {
    if (o.candidates < int(h.nvec) || o.candidates > 64 || o.iterations < 1 ||
        o.restart < 2 || o.restart > 24 || o.warm < 0 || o.batch < 0 ||
        (o.batch && o.batch > 32) || !(o.tolerance > 0 && o.tolerance < 1))
      throw std::runtime_error("invalid near-null options");
    size_t n = size_t(h.rows) * h.inbs;
    uint32_t carry = h.nvec + 1;
    uint32_t stride = o.batch
                          ? std::min(uint32_t(o.candidates), carry + o.batch)
                          : uint32_t(o.candidates);
    uint32_t represented = 0;
    std::vector<Cx<float>> candidates(n * stride);
    auto started = Clock::now();
    vectors_log << "candidate,iterations,relative_Sv_to_Sy,Sv_over_v,seconds\n";
    {
      struct Diagonal {
        Schur &op;
        void encode(id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                    id<MTLBuffer> y) {
          op.block(e, op.Aodd, x, y);
        }
        void check_coarse() {}
      } diagonal{op};
      GPUGenericFGMRES<Schur, Diagonal> solver(uint32_t(n), 1, op, diagonal,
                                               o.restart);
      solver.allow_inexact = true;
      auto alloc = [&]() {
        auto b = [op.m.device newBufferWithLength:n * 8
                                          options:MTLResourceStorageModeShared];
        if (!b)
          throw std::runtime_error("near-null vector allocation");
        return b;
      };
      auto y = alloc(), sy = alloc();
      TuneVector update(op.m);
      auto norm = [&](id<MTLBuffer> b) {
        auto *p = static_cast<const Cx<float> *>([b contents]);
        std::vector<double> sums(h.workers);
        chunks(uint32_t(n), h.workers,
               [&](int worker, uint32_t first, uint32_t last) {
                 double sum = 0;
                 for (size_t i = first; i < last; ++i)
                   sum += double(p[i].re) * p[i].re + double(p[i].im) * p[i].im;
                 sums[worker] = sum;
               });
        double sum = 0;
        for (double v : sums)
          sum += v;
        return std::sqrt(sum);
      };
      for (int v = 0; v < o.candidates; ++v) {
        auto begin = Clock::now();
        auto *yp = static_cast<Cx<float> *>([y contents]);
        double scale = 1 / std::sqrt(double(n));
        chunks(uint32_t(n), h.workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = first; i < last; ++i) {
            auto z = chroma_bridge::source(i, 771337 + v);
            yp[i] = {float(scale * z.re), float(scale * z.im)};
          }
        });
        for (int first = 0; first < o.warm; first += 8) {
          @autoreleasepool {
            auto cmd = [op.m.queue commandBuffer];
            auto e = [cmd computeCommandEncoder];
            for (int it = first; it < std::min(o.warm, first + 8); ++it) {
              op.encode(e, y, sy);
              op.block(e, op.Aodd, sy, op.inverted);
              update.encode(e, op.inverted, y, uint32_t(n), -.65f, 1);
            }
            solver.finish(cmd, e);
          }
        }
        {
          auto cmd = [op.m.queue commandBuffer];
          auto e = [cmd computeCommandEncoder];
          op.encode(e, y, sy);
          solver.finish(cmd, e);
        }
        double initial = norm(sy);
        auto info = solver.run_device(sy, o.tolerance, o.iterations);
        {
          auto cmd = [op.m.queue commandBuffer];
          auto e = [cmd computeCommandEncoder];
          update.encode(e, solver.delta, y, uint32_t(n), -1, 1);
          op.encode(e, y, sy);
          solver.finish(cmd, e);
        }
        double size = norm(y), product = norm(sy);
        if (!(size > 0) || !std::isfinite(product) || product > initial * 1.01)
          throw std::runtime_error(
              "homogeneous near-null solve failed to reduce S v");
        chunks(uint32_t(n), h.workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = first; i < last; ++i)
            candidates[i * stride + represented] = {float(yp[i].re / size),
                                                    float(yp[i].im / size)};
        });
        ++represented;
        if (o.batch && represented == stride && v + 1 < o.candidates) {
          auto compressed = compress_near_null_svd(h, candidates, represented,
                                                   nullptr, stride, carry);
          represented = carry;
          std::cout << "NEAR_NULL_CARRY generated=" << v + 1
                    << " retained=" << carry
                    << " represented_energy=" << compressed.kept_energy
                    << " fraction=" << compressed.retained_energy << std::endl;
        }
        vectors_log << v << ',' << info.max_iterations << ','
                    << product / initial << ',' << product / size << ','
                    << seconds(begin) << '\n';
        vectors_log.flush();
        std::cout << "NEAR_NULL candidate=" << v << " K=" << info.max_iterations
                  << " residual_reduction=" << product / initial
                  << " Sv_over_v=" << product / size
                  << " seconds=" << seconds(begin) << std::endl;
      }
    }
    auto before_svd = Clock::now();
    auto stats =
        compress_near_null_svd(h, candidates, represented, &spectrum, stride);
    std::cout << "NEAR_NULL_SVD seconds=" << seconds(started)
              << " compression_seconds=" << seconds(before_svd)
              << " candidates=" << o.candidates << " retained=" << h.nvec
              << " batch=" << o.batch << " pool_columns=" << stride
              << " energy_lower_bound=" << stats.kept_energy / o.candidates
              << " energy=" << stats.retained_energy
              << " minimum_energy=" << stats.minimum_retained_energy
              << " minimum_rank_ratio=" << stats.minimum_rank_ratio
              << " orthogonality=" << stats.maximum_orthogonality_error
              << std::endl;
  }

} // namespace MetalCloverEnv
