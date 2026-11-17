// Self-contained public-API example: identity links and manufactured full RHS.
// No Chroma installation, gauge files, or internal implementation headers
// needed.
#include <cmath>
#include <cstdio>
#include <memory>
#include <metalclover/metalclover.h>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace MetalCloverEnv;

namespace {
  void check(int status) {
    if (status)
      throw std::runtime_error(metalclover_last_error());
  }

  size_t count_arg(const char *text) {
    std::string value(text);
    if (value.empty() ||
        value.find_first_not_of("0123456789") != std::string::npos)
      throw std::runtime_error(
          "RHS count and batch cap must be integers in 1..64");
    size_t result = std::stoul(value);
    if (result < 1 || result > 64)
      throw std::runtime_error("RHS count and batch cap must be in 1..64");
    return result;
  }

  double relative_error(const std::vector<double> &a,
                        const std::vector<double> &b) {
    double error = 0, scale = 0;
    for (size_t i = 0; i < a.size(); ++i) {
      error += (a[i] - b[i]) * (a[i] - b[i]);
      scale += b[i] * b[i];
    }
    return std::sqrt(error / (scale > 0 ? scale : 1));
  }
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::printf("Usage: %s [bicg|mg] [rhs_count=5] [batch_cap=3]\n"
                  "Counts are 1..64. Uses a synthetic 4^4 identity gauge.\n",
                  argv[0]);
      return 0;
    }
    const std::string method = argc > 1 ? argv[1] : "bicg";
    if (argc > 4 || (method != "bicg" && method != "mg"))
      throw std::runtime_error(
          "Usage: multi_rhs [bicg|mg] [rhs_count] [batch_cap]");
    const size_t rhs = argc > 2 ? count_arg(argv[2]) : 5;
    const size_t cap = argc > 3 ? count_arg(argv[3]) : 3;

    metalclover_options options;
    metalclover_default_options(&options);
    size_t volume = 1;
    for (int mu = 0; mu < 4; ++mu) {
      options.shape[mu] = 4;
      options.blocks[mu] = 2; // First coarse extents must be even for MG.
      volume *= options.shape[mu];
    }
    options.mass = 0.5;
    options.method = method == "mg" ? METALCLOVER_MG : METALCLOVER_BICGSTAB;
    // On this tiny free field, long homogeneous filters can collapse the
    // sample rank. Keep the second-level filter short for this demonstration.
    options.second_filter = 1;

    // Public layout: interleaved doubles; links are site,mu,row,col.
    // Identity links need no smearing. The default temporal AP sign is applied
    // by the library, and the default action coefficients are isotropic.
    std::vector<double> gauge(72 * volume, 0);
    for (size_t site = 0; site < volume; ++site)
      for (size_t mu = 0; mu < 4; ++mu)
        for (size_t color = 0; color < 3; ++color)
          gauge[2 * (((site * 4 + mu) * 3 + color) * 3 + color)] = 1;
    metalclover_context *raw = nullptr;
    check(metalclover_create(&options, gauge.data(), &raw));
    std::unique_ptr<metalclover_context, decltype(&metalclover_destroy)>
        context(raw, metalclover_destroy);
    check(metalclover_set_rhs_batch(context.get(), cap));

    const size_t doubles = 24 * volume;
    std::vector<std::vector<double>> exact(rhs, std::vector<double>(doubles));
    std::vector<std::vector<double>> sources(rhs, std::vector<double>(doubles));
    std::vector<std::vector<double>> solutions(rhs,
                                               std::vector<double>(doubles, 0));
    std::vector<const void *> b(rhs);
    std::vector<void *> x(rhs);
    for (size_t r = 0; r < rhs; ++r) {
      std::mt19937 random(static_cast<unsigned>(771337 + r));
      for (double &value : exact[r])
        value = double(random()) / double(std::mt19937::max()) - 0.5;
      // Manufacture b=M*x_exact using the public CPU FP64 operator helper.
      check(metalclover_apply_full(context.get(), exact[r].data(),
                                   sources[r].data()));
      b[r] = sources[r].data();
      x[r] = solutions[r].data();
    }

    std::vector<metalclover_solve_info> info(rhs);
    size_t failed = SIZE_MAX;
    const double tolerance = 1e-12;
    const int status = metalclover_solve_full_array(
        context.get(), rhs, b.data(), x.data(), tolerance, 2000, nullptr,
        nullptr, info.data(), &failed);
    if (status) {
      // Save this before another API call clears the thread-local error.
      const std::string error = metalclover_last_error();
      throw std::runtime_error("array failed at RHS " + std::to_string(failed) +
                               ": " + error);
    }
    metalclover_batch_info batch{};
    check(metalclover_get_batch_info(context.get(), &batch));
    std::printf("method=%s RHS=%zu cap=%zu max_resident=%zu tiles=%zu "
                "setup=%.6fs solve=%.6fs\n",
                method.c_str(), rhs, cap, batch.maximum_resident_rhs,
                batch.tiles, metalclover_setup_seconds(context.get()),
                batch.seconds);

    // Recompute full residuals and compare with the known solutions. This is
    // a public-API consistency check, not an independent Chroma validation.
    std::vector<double> product(doubles);
    for (size_t r = 0; r < rhs; ++r) {
      check(metalclover_apply_full(context.get(), solutions[r].data(),
                                   product.data()));
      const double residual = relative_error(product, sources[r]);
      const double error = relative_error(solutions[r], exact[r]);
      std::printf("RHS %zu: converged=%d iterations=%d residual=%.3e "
                  "solution_error=%.3e\n",
                  r, info[r].converged, info[r].krylov_iterations, residual,
                  error);
      if (!info[r].converged || !std::isfinite(residual) ||
          !std::isfinite(error) || residual > tolerance || error > 1e-9)
        throw std::runtime_error(
            "full residual or known-solution check failed");
    }
    std::puts("PASS: all multi-RHS full-system checks passed");
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "multi_rhs: %s\n", error.what());
    return 1;
  }
}
