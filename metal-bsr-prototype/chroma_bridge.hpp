#pragma once

#include <array>
#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace MetalCloverEnv {

  // Single-node interchange. All scalars are little-endian on Apple Silicon.
  // Link order: x-fastest site, direction, row, column, real/imaginary.
  // Spinor order: x-fastest site, color, spin, real/imaginary.
  namespace chroma_bridge {
    inline constexpr std::array<char, 8> magic{'C', 'H', 'R', 'M',
                                               'G', '0', '0', '1'};
    inline constexpr uint32_t probe_rhs = 0xffffffffu;

    struct Complex {
      double re, im;
    };
    static_assert(sizeof(Complex) == 2 * sizeof(double));

    inline uint32_t volume(uint32_t extent) {
      if (extent < 4 || extent > 32 || extent % 2)
        throw std::runtime_error(
            "Chroma bridge requires an even hypercube with extent 4–32");
      return extent * extent * extent * extent;
    }

    inline Complex source(uint64_t component, uint32_t rhs) {
      uint64_t z = component + (uint64_t(rhs) << 32) + 0x9e3779b97f4a7c15ULL;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      switch ((z ^ (z >> 31)) & 3) {
      case 0:
        return {1, 0};
      case 1:
        return {-1, 0};
      case 2:
        return {0, 1};
      default:
        return {0, -1};
      }
    }

    template <class T>
    void write(std::ostream &out, const T *data, size_t count) {
      if (!out.write(reinterpret_cast<const char *>(data), count * sizeof(T)))
        throw std::runtime_error("failed to write Chroma bridge payload");
    }

    template <class T> void read(std::istream &in, T *data, size_t count) {
      if (!in.read(reinterpret_cast<char *>(data), count * sizeof(T)))
        throw std::runtime_error("truncated Chroma bridge payload");
    }
  } // namespace chroma_bridge

} // namespace MetalCloverEnv
