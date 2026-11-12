#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace MetalCloverEnv {

  using LatticeShape = std::array<uint32_t, 4>;

  inline LatticeShape hypercube(uint32_t extent) {
    return {extent, extent, extent, extent};
  }

  inline uint32_t lattice_volume(const LatticeShape &shape) {
    uint64_t volume = 1;
    for (uint32_t n : shape) {
      if (n < 2 || n % 2)
        throw std::runtime_error(
            "lattice extents must be positive even integers");
      volume *= n;
    }
    if (volume > std::numeric_limits<uint32_t>::max() / 12)
      throw std::runtime_error("lattice exceeds component index range");
    return uint32_t(volume);
  }

} // namespace MetalCloverEnv
