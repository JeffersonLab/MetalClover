#include "chroma_import.hpp"
#include "compact_residual.hpp"
#include <algorithm>
#include <iomanip>
#include <unistd.h>

struct ResidualTimings {
    double operator_ms = 0;
    double subtract_norm_ms = 0;
};

double compare_residual(const std::vector<D>& reference, const std::vector<D>& result,
                        const std::vector<D>& rhs) {
    long double error = 0, norm = 0;
    for (size_t i = 0; i < reference.size(); ++i) {
        auto delta = minus_z(reference[i], result[i]);
        error += (long double)delta.re * delta.re + (long double)delta.im * delta.im;
        norm += (long double)rhs[i].re * rhs[i].re + (long double)rhs[i].im * rhs[i].im;
    }
    double relative = std::sqrt(error / std::max(norm, 1e-300L));
    if (!std::isfinite(relative) || relative > 5e-13)
        throw std::runtime_error("optimized FP64 residual differs from original stencil");
    return relative;
}

// Original full residual path, retained as the numerical and timing reference.
ResidualTimings reference_residual(const CompactClover& c, const std::vector<D>& x,
                                  std::vector<D>& residual, std::vector<double>& sums,
                                  int workers) {
    const auto& a = c.w.a;
    auto start = Clock::now();
    apply_cpu(a, x.data(), residual.data(), workers);
    double operator_ms = seconds(start) * 1000;
    start = Clock::now();
    std::fill(sums.begin(), sums.end(), 0.);
    chunks(a.p.rows * 12, workers, [&](int worker, uint32_t first, uint32_t last) {
        auto* partial = sums.data() + size_t(worker) * a.p.rhs;
        for (size_t q = first; q < last; ++q)
            for (uint32_t r = 0; r < a.p.rhs; ++r) {
                size_t j = q * a.p.rhs + r;
                residual[j] = minus_z(a.b[j], residual[j]);
                auto z = residual[j];
                partial[r] += z.re * z.re + z.im * z.im;
            }
    });
    return {operator_ms, seconds(start) * 1000};
}

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3 && argc != 5)
            throw std::runtime_error("Usage: residual_profile BRIDGE [ITERATIONS [SOLUTIONS RHS]]");
        int iterations = argc >= 3 ? std::stoi(argv[2]) : 7;
        if (iterations < 1 || iterations > 1000) throw std::runtime_error("invalid iterations");
        ChromaImport source(argv[1]);
        CompactClover c(source.extent, source.mass, 8, .25, source.links());
        source.validate(c, 8);
        c.rhs(4, false);
        c.w.a.b.resize(c.w.a.p.outputs());
        CompactOuter outer(c.w.a, 8);
        std::ifstream no_file;
        set_chroma_rhs(c, outer, 0, 8, no_file);
        std::vector<D> x(outer.x.size()), residual(x.size()), optimized(x.size());
        for (size_t j = 0; j < x.size(); ++j) x[j] = compact_truth(j / 4, j % 4);
        std::cout << "PROFILE pid=" << getpid() << " L=" << source.extent
                  << " rhs=4 iterations=" << iterations << std::endl;
        reference_residual(c, x, residual, outer.sums, 8);
        std::vector<double> operators, reductions, fused;
        double error = 0;
        auto optimized_run = [&] {
            auto start = Clock::now();
            compact_fp64::residual(c, x.data(), c.w.a.b.data(), optimized.data(),
                                   outer.sums.data(), 8);
            fused.push_back(seconds(start) * 1000);
        };
        optimized_run();
        fused.clear();
        for (int i = 0; i < iterations; ++i) {
            if (i % 2 == 0) optimized_run();
            auto t = reference_residual(c, x, residual, outer.sums, 8);
            if (i % 2 != 0) optimized_run();
            operators.push_back(t.operator_ms);
            reductions.push_back(t.subtract_norm_ms);
            error = std::max(error, compare_residual(residual, optimized, c.w.a.b));
        }
        auto median = [](std::vector<double> a) {
            std::sort(a.begin(), a.end());
            return a[a.size() / 2];
        };
        std::cout << std::setprecision(10)
                  << "PROFILE baseline_operator_median_ms=" << median(operators)
                  << " baseline_subtract_norm_median_ms=" << median(reductions)
                  << " optimized_median_ms=" << median(fused)
                  << " residual_difference_over_b=" << error
                  << " checksum=" << outer.sums[0] << std::endl;

        if (argc == 5) {
            std::ifstream solutions(argv[3], std::ios::binary), rhs_file(argv[4], std::ios::binary);
            std::vector<chroma_bridge::Complex> buffer(x.size() / 4);
            for (uint32_t r = 0; r < 4; ++r) {
                chroma_bridge::read(solutions, buffer.data(), buffer.size());
                for (size_t q = 0; q < buffer.size(); ++q)
                    x[q * 4 + r] = {buffer[q].re, buffer[q].im};
                chroma_bridge::read(rhs_file, buffer.data(), buffer.size());
                for (size_t q = 0; q < buffer.size(); ++q)
                    c.w.a.b[q * 4 + r] = {buffer[q].re, buffer[q].im};
            }
            compact_rhs_norms(c.w.a, outer, 8);
            reference_residual(c, x, residual, outer.sums, 8);
            compact_fp64::residual(c, x.data(), c.w.a.b.data(), optimized.data(),
                                   outer.sums.data(), 8);
            double difference = compare_residual(residual, optimized, c.w.a.b);
            double worst = 0;
            for (uint32_t r = 0; r < 4; ++r) {
                double sum = 0;
                for (int worker = 0; worker < 8; ++worker) sum += outer.sums[worker * 4 + r];
                worst = std::max(worst, std::sqrt(sum) / c.w.a.bnorm[r]);
            }
            if (!std::isfinite(worst) || worst > 1e-10)
                throw std::runtime_error("optimized check rejects saved Chroma solution");
            std::cout << "CANCELLATION PASS residual_difference_over_b=" << difference
                      << " worst_full_relative=" << worst << std::endl;
        }

        if (source.extent == 4) {
            for (uint32_t rhs : {1, 2, 3, 5, 8, 12}) {
                c.rhs(rhs, false);
                c.w.a.b.assign(c.w.a.p.outputs(), {1, -.5});
                std::vector<D> input(c.w.a.b.size()), ref(input.size()), got(input.size());
                for (size_t q = 0; q < input.size(); ++q) input[q] = compact_truth(q / rhs, q % rhs);
                for (int workers : {1, 8}) {
                    std::vector<double> sums(size_t(workers) * rhs);
                    reference_residual(c, input, ref, sums, workers);
                    auto ref_sums = sums;
                    compact_fp64::residual(c, input.data(), c.w.a.b.data(), got.data(),
                                           sums.data(), workers);
                    compare_residual(ref, got, c.w.a.b);
                    for (uint32_t r = 0; r < rhs; ++r) {
                        double old_norm = 0, new_norm = 0;
                        for (int worker = 0; worker < workers; ++worker) {
                            old_norm += ref_sums[worker * rhs + r];
                            new_norm += sums[worker * rhs + r];
                        }
                        if (std::abs(old_norm - new_norm) > 5e-13 * old_norm)
                            throw std::runtime_error("optimized per-RHS norm mismatch");
                    }
                }
            }
            std::cout << "PASS RHS tails and one/eight-worker norms" << std::endl;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
