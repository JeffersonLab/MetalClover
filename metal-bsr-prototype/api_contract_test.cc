#include "include/metalclover/metalclover.h"
#include "production_bridge.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace MetalCloverEnv;
namespace {
  using Z = chroma_bridge::Complex;
  void check(int code) {
    if (code)
      throw std::runtime_error(metalclover_last_error());
  }
  double relative(const std::vector<Z> &a, const std::vector<Z> &b) {
    double e = 0, n = 0;
    for (size_t i = 0; i < a.size(); ++i) {
      double re = a[i].re - b[i].re, im = a[i].im - b[i].im;
      e += re * re + im * im;
      n += b[i].re * b[i].re + b[i].im * b[i].im;
    }
    return std::sqrt(e / (n > 0 ? n : 1));
  }
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 2 && argc != 3)
      throw std::runtime_error(
          "Usage: api_contract_test ALIGNED_BRIDGE [bicg|mg]");
    std::ifstream in(argv[1], std::ios::binary);
    auto h = production::read_header(in);
    size_t n = size_t(h.volume()) * 12;
    std::vector<Z> links(size_t(h.volume()) * 36), b(n), x(n), ax(n),
        odd(n / 2), reduced(n / 2), reconstructed(n);
    chroma_bridge::read(in, links.data(), links.size());
    metalclover_options o;
    metalclover_default_options(&o);
    std::copy(h.shape.begin(), h.shape.end(), o.shape);
    o.mass = h.mass;
    o.xi0 = h.xi0;
    o.nu = h.nu;
    o.clover_r = h.clover_r;
    o.clover_t = h.clover_t;
    o.method = argc == 3 && std::string(argv[2]) == "mg" ? METALCLOVER_MG
                                                         : METALCLOVER_BICGSTAB;
    for (int d = 0; d < 4; ++d)
      if (h.shape[d] % o.blocks[d] || (h.shape[d] / o.blocks[d]) % 2)
        o.blocks[d] = h.shape[d] / 4;
    o.second_filter = 16;
    o.bottom_iterations = 4;
    metalclover_context *raw = nullptr;
    check(metalclover_create(&o, links.data(), &raw));
    std::unique_ptr<metalclover_context, decltype(&metalclover_destroy)>
        context(raw, metalclover_destroy);
    check(metalclover_set_rhs_batch(raw, 3));
    if (!metalclover_set_rhs_batch(raw, 0) ||
        !metalclover_set_rhs_batch(raw, 65))
      throw std::runtime_error("invalid batch capacity accepted");
    for (size_t i = 0; i < n; ++i)
      b[i] = chroma_bridge::source(i, 20261112);
    check(metalclover_reduce_rhs(raw, b.data(), reduced.data()));
    metalclover_solve_info full_info{}, odd_info{};
    check(metalclover_solve_full(raw, b.data(), x.data(), 1e-12, 65536, nullptr,
                                 nullptr, &full_info));
    check(metalclover_solve_odd(raw, reduced.data(), odd.data(), 1e-12, 65536,
                                nullptr, nullptr, &odd_info));
    check(metalclover_reconstruct(raw, b.data(), odd.data(),
                                  reconstructed.data()));
    double agreement = relative(reconstructed, x);
    if (agreement > 1e-9)
      throw std::runtime_error("full/odd API solutions differ");
    check(metalclover_apply_full(raw, reconstructed.data(), ax.data()));
    double residual = relative(ax, b);
    if (residual > 1e-12)
      throw std::runtime_error("reconstructed full API residual failed");
    metalclover_solve_info exact{};
    check(metalclover_solve_full(raw, b.data(), x.data(), 1e-12, 65536, nullptr,
                                 nullptr, &exact));
    if (exact.krylov_iterations || exact.nominal_operator_flops)
      throw std::runtime_error("exact initial guess contract failed");
    std::fill(x.begin(), x.end(), Z{});
    metalclover_solve_info failure{};
    if (!metalclover_solve_full(raw, b.data(), x.data(), 1e-12, 1, nullptr,
                                nullptr, &failure))
      throw std::runtime_error("one-iteration budget incorrectly accepted");
    if (metalclover_solve_full(raw, b.data(), x.data(), 0, 65536, nullptr,
                               nullptr, &failure) == 0)
      throw std::runtime_error("invalid tolerance accepted");
    if (!(full_info.fine_schur_applications > 0) ||
        !(odd_info.nominal_operator_flops > 0) ||
        std::abs(full_info.nominal_operator_gflops * full_info.seconds * 1e9 /
                     full_info.nominal_operator_flops -
                 1) > 1e-12)
      throw std::runtime_error("solver info accounting failed");
    {
      std::vector<Z> b1(n), full0(n), full1(n), odd0(n / 2), odd1(n / 2),
          reduced1(n / 2), rebuilt(n);
      for (size_t i = 0; i < n; ++i)
        b1[i] = chroma_bridge::source(i, 20261113);
      check(metalclover_reduce_rhs(raw, b1.data(), reduced1.data()));
      const void *sources[] = {b.data(), b1.data()};
      void *solutions[] = {full0.data(), full1.data()};
      metalclover_solve_info batch[2]{};
      size_t failed = 0;
      check(metalclover_solve_full_array(raw, 2, sources, solutions, 1e-12,
                                         65536, nullptr, nullptr, batch,
                                         &failed));
      if (failed != SIZE_MAX || !batch[0].converged || !batch[1].converged)
        throw std::runtime_error("full array result contract");
      const void *odd_sources[] = {reduced.data(), reduced1.data()};
      void *odd_solutions[] = {odd0.data(), odd1.data()};
      check(metalclover_solve_odd_array(raw, 2, odd_sources, odd_solutions,
                                        1e-12, 65536, nullptr, nullptr, batch,
                                        &failed));
      check(
          metalclover_reconstruct(raw, b1.data(), odd1.data(), rebuilt.data()));
      if (relative(rebuilt, full1) > 1e-9 ||
          relative(full0, reconstructed) > 1e-9)
        throw std::runtime_error("array/scalar/full/odd solution disagreement");
      check(metalclover_solve_full_array(raw, 0, nullptr, nullptr, 1e-12, 65536,
                                         nullptr, nullptr, nullptr, &failed));
      if (failed != SIZE_MAX)
        throw std::runtime_error("empty array result contract");
      void *aliases[] = {full0.data(), full0.data()};
      if (!metalclover_solve_full_array(raw, 2, sources, aliases, 1e-12, 65536,
                                        nullptr, nullptr, batch, &failed))
        throw std::runtime_error("aliased array outputs accepted");
      // Unequal scales, zero/exact guesses, collinear sources, a nonzero guess,
      // and a partial final tile exercise independent lane state and FP64
      // checks.
      constexpr size_t count = 5;
      std::vector<std::vector<Z>> many_b(count, std::vector<Z>(n)),
          many_x(count, std::vector<Z>(n)), odd_b(count, std::vector<Z>(n / 2)),
          odd_x(count, std::vector<Z>(n / 2));
      const double scale[count] = {1e-8, 0, 1e8, 1, 1};
      for (size_t col = 0; col < count; ++col)
        for (size_t q = 0; q < n; ++q) {
          auto z = col == 4 ? b1[q] : b[q];
          many_b[col][q] = {scale[col] * z.re, scale[col] * z.im};
          if (col == 3)
            many_x[col][q] = full0[q];
          if (col == 4)
            many_x[col][q] = {.01 * b[q].re, .01 * b[q].im};
        }
      std::vector<const void *> in_many(count), in_odd(count);
      std::vector<void *> out_many(count), out_odd(count);
      std::vector<metalclover_solve_info> info(count);
      for (size_t col = 0; col < count; ++col) {
        in_many[col] = many_b[col].data();
        out_many[col] = many_x[col].data();
        check(
            metalclover_reduce_rhs(raw, many_b[col].data(), odd_b[col].data()));
        in_odd[col] = odd_b[col].data();
        out_odd[col] = odd_x[col].data();
      }
      check(metalclover_solve_full_array(raw, count, in_many.data(),
                                         out_many.data(), 1e-12, 65536, nullptr,
                                         nullptr, info.data(), &failed));
      metalclover_batch_info summary{};
      check(metalclover_get_batch_info(raw, &summary));
      if (summary.maximum_resident_rhs != 3 || summary.tiles != 2 ||
          summary.rhs_count != count)
        throw std::runtime_error("true GPU tiling metadata failed");
      uint64_t work = 0;
      double elapsed = 0;
      for (size_t col = 0; col < count; ++col) {
        check(metalclover_apply_full(raw, many_x[col].data(), ax.data()));
        if (!(relative(ax, many_b[col]) <= 1e-12) || !info[col].converged)
          throw std::runtime_error(
              "mixed-scale batch full FP64 residual failed");
        work += info[col].nominal_operator_flops;
        elapsed += info[col].seconds;
        if (col == 1 || col == 3) {
          if (info[col].krylov_iterations || info[col].nominal_operator_flops)
            throw std::runtime_error("zero/exact batch lane did work");
        }
        if (col != 1) {
          std::vector<Z> expected(n);
          for (size_t q = 0; q < n; ++q) {
            auto z = col == 4 ? full1[q] : full0[q];
            expected[q] = {scale[col] * z.re, scale[col] * z.im};
          }
          if (relative(many_x[col], expected) > 1e-8)
            throw std::runtime_error(
                "batch/scalar mixed-scale solutions disagree");
        }
      }
      if (work != summary.nominal_operator_flops ||
          std::abs(elapsed / summary.seconds - 1) > 1e-12)
        throw std::runtime_error("aggregate batch accounting failed");
      check(metalclover_solve_odd_array(raw, count, in_odd.data(),
                                        out_odd.data(), 1e-12, 65536, nullptr,
                                        nullptr, info.data(), &failed));
      for (size_t col = 0; col < count; ++col) {
        check(metalclover_reconstruct(raw, many_b[col].data(),
                                      odd_x[col].data(), rebuilt.data()));
        check(metalclover_apply_full(raw, rebuilt.data(), ax.data()));
        if (!(relative(ax, many_b[col]) <= 1e-12) ||
            relative(rebuilt, many_x[col]) > 1e-8)
          throw std::runtime_error("batch odd/full reconstruction disagrees");
      }
      for (auto &v : many_x)
        std::fill(v.begin(), v.end(), Z{});
      if (!metalclover_solve_full_array(raw, count, in_many.data(),
                                        out_many.data(), 1e-12, 1, nullptr,
                                        nullptr, info.data(), &failed))
        throw std::runtime_error(
            "batch total iteration budget accepted incorrectly");
      for (const auto &v : info)
        if (v.krylov_iterations > 1)
          throw std::runtime_error(
              "per-column batch iteration budget exceeded");
      check(metalclover_set_rhs_batch(raw, 1));
      std::fill(full0.begin(), full0.end(), Z{});
      std::fill(full1.begin(), full1.end(), Z{});
      check(metalclover_solve_full_array(raw, 2, sources, solutions, 1e-12,
                                         65536, nullptr, nullptr, batch,
                                         &failed));
      check(metalclover_get_batch_info(raw, &summary));
      if (summary.maximum_resident_rhs != 1 || summary.tiles != 2)
        throw std::runtime_error("scalar batch-cap fallback failed");
      void *overlap[] = {full0.data(), full0.data() + 1};
      if (!metalclover_solve_full_array(raw, 2, sources, overlap, 1e-12, 65536,
                                        nullptr, nullptr, batch, &failed))
        throw std::runtime_error("overlapping field ranges accepted");
      std::cout << "PASS mixed-scale zero/exact/collinear/nonzero-guess GPU "
                   "tiles=3+2, "
                   "full/odd FP64 checks, scalar agreement, per-column budgets "
                   "and batch accounting"
                << std::endl;
    }
    // Invalid context construction must not publish a partially built handle.
    metalclover_context *invalid = reinterpret_cast<metalclover_context *>(1);
    auto bad = o;
    bad.abi_version = 999;
    if (!metalclover_create(&bad, links.data(), &invalid) || invalid)
      throw std::runtime_error("invalid ABI context accepted");
    bad = o;
    bad.mass = -10;
    invalid = nullptr;
    if (!metalclover_create(&bad, links.data(), &invalid) || invalid)
      throw std::runtime_error("nonpositive local LDL pivot accepted");
    std::cout << "PASS C ABI full/odd/reconstruction, exact guess, "
                 "budget/tolerance/ABI errors, "
                 "structured work info; reconstructed_relative="
              << residual << " solution_agreement=" << agreement << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR " << e.what() << std::endl;
    return 1;
  }
}
