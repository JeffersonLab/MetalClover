#pragma once
#include "chroma_bridge.hpp"
#include "compact_solve.hpp"

namespace MetalCloverEnv {

  class ChromaImport {
    std::ifstream input;

    double compare(const std::vector<D> &values) {
      constexpr size_t block = 4096;
      std::array<chroma_bridge::Complex, block> expected;
      double error = 0, norm = 0;
      for (size_t first = 0; first < values.size(); first += block) {
        size_t n = std::min(block, values.size() - first);
        chroma_bridge::read(input, expected.data(), n);
        for (size_t j = 0; j < n; ++j) {
          auto ref = expected[j];
          auto got = values[first + j];
          error += std::norm(CD(got.re - ref.re, got.im - ref.im));
          norm += std::norm(CD(ref.re, ref.im));
        }
      }
      return std::sqrt(error / norm);
    }

  public:
    uint32_t extent;
    double mass;

    explicit ChromaImport(const std::string &path)
        : input(path, std::ios::binary) {
      if (!input)
        throw std::runtime_error("cannot open Chroma bridge file");
      std::array<char, chroma_bridge::magic.size()> magic;
      chroma_bridge::read(input, magic.data(), magic.size());
      if (magic != chroma_bridge::magic)
        throw std::runtime_error("unsupported Chroma bridge format");
      chroma_bridge::read(input, &extent, 1);
      chroma_bridge::read(input, &mass, 1);
      chroma_bridge::volume(extent);
      if (!std::isfinite(mass))
        throw std::runtime_error("nonfinite Chroma mass");
    }

    std::vector<M3> links() {
      std::vector<M3> gauge(size_t(chroma_bridge::volume(extent)) * 4);
      constexpr size_t block = 256;
      std::array<chroma_bridge::Complex, 9 * block> raw;
      for (size_t first = 0; first < gauge.size(); first += block) {
        size_t n = std::min(block, gauge.size() - first);
        chroma_bridge::read(input, raw.data(), 9 * n);
        for (size_t i = 0; i < n; ++i)
          for (size_t j = 0; j < 9; ++j) {
            auto z = raw[i * 9 + j];
            if (!std::isfinite(z.re) || !std::isfinite(z.im))
              throw std::runtime_error("nonfinite Chroma gauge link");
            gauge[first + i][j] = {z.re, z.im};
          }
      }
      return gauge;
    }

    void validate(CompactClover &c, int workers) {
      c.rhs(1, false);
      std::vector<D> x(c.w.a.p.outputs()), result(x.size());
      for (size_t q = 0; q < x.size(); ++q) {
        auto z = chroma_bridge::source(q, chroma_bridge::probe_rhs);
        x[q] = {z.re, z.im};
      }
      apply_cpu(c.w.a, x.data(), result.data(), workers);
      double full = compare(result);
      result.resize(c.parity.oe.p.outputs());
      for (size_t i = 0; i < c.parity.odd.size(); ++i)
        std::copy_n(x.data() + size_t(c.parity.odd[i]) * 12, 12,
                    x.data() + i * 12);
      CPUCloverSchur<double, CompactClover> schur(c, workers);
      schur.apply(x.data(), result.data());
      double odd = compare(result);
      if (input.peek() != EOF)
        throw std::runtime_error("unexpected bytes after Chroma references");
      std::cout << "CHROMA FP64 full_relative=" << full
                << " odd_schur_relative=" << odd << std::endl;
      if (!std::isfinite(full) || !std::isfinite(odd) || full > 1e-12 ||
          odd > 1e-12)
        throw std::runtime_error("standalone operator disagrees with Chroma");
    }
  };

  inline void set_chroma_rhs(CompactClover &c, CompactOuter &outer,
                             uint32_t offset, int workers,
                             std::ifstream &source_file) {
    auto &a = c.w.a;
    uint32_t rhs = a.p.rhs;
    const size_t components = size_t(a.p.rows) * 12;
    if (source_file.is_open()) {
      constexpr size_t block = 4096;
      std::array<chroma_bridge::Complex, block> buffer;
      source_file.seekg(std::streamoff(offset) * components *
                        sizeof(buffer[0]));
      for (uint32_t r = 0; r < rhs; ++r)
        for (size_t first = 0; first < components; first += block) {
          size_t n = std::min(block, components - first);
          chroma_bridge::read(source_file, buffer.data(), n);
          for (size_t j = 0; j < n; ++j) {
            auto z = buffer[j];
            if (!std::isfinite(z.re) || !std::isfinite(z.im))
              throw std::runtime_error("nonfinite Chroma RHS");
            a.b[(first + j) * rhs + r] = {z.re, z.im};
          }
        }
    } else {
      chunks(components, workers, [&](int, uint32_t first, uint32_t last) {
        for (size_t q = first; q < last; ++q)
          for (uint32_t r = 0; r < rhs; ++r) {
            auto z = chroma_bridge::source(q, offset + r);
            a.b[q * rhs + r] = {z.re, z.im};
          }
      });
    }
    compact_rhs_norms(a, outer, workers);
    for (double norm : a.bnorm)
      if (!std::isfinite(norm))
        throw std::runtime_error("nonfinite Chroma RHS norm");
  }

  inline void write_chroma_solutions(std::ostream &out,
                                     const CompactOuter &outer, uint32_t rhs) {
    constexpr size_t block = 4096;
    std::array<chroma_bridge::Complex, block> buffer;
    const size_t components = outer.x.size() / rhs;
    for (uint32_t r = 0; r < rhs; ++r)
      for (size_t first = 0; first < components; first += block) {
        size_t n = std::min(block, components - first);
        for (size_t j = 0; j < n; ++j) {
          auto z = outer.x[(first + j) * rhs + r];
          buffer[j] = {z.re, z.im};
        }
        chroma_bridge::write(out, buffer.data(), n);
      }
    out.flush();
    if (!out)
      throw std::runtime_error("failed to flush Chroma solutions");
  }

} // namespace MetalCloverEnv
