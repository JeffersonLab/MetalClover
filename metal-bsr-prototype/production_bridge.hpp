#pragma once
#include "chroma_bridge.hpp"
#include "lattice_geometry.hpp"
#include <cmath>

namespace production {
  inline constexpr std::array<char, 8> magic{'C', 'H', 'R', 'M',
                                             'G', '0', '0', '2'};
  inline constexpr std::array<char, 8> eigen_magic{'C', 'H', 'R', 'E',
                                                   'V', '0', '0', '1'};
  struct Header {
    LatticeShape shape{};
    double mass = 0, xi0 = 1, nu = 1, clover_r = 1, clover_t = 1;
    uint32_t source_time = 0, first_time = 0, time_count = 0;
    uint32_t volume() const { return lattice_volume(shape); }
    void validate() const {
      volume();
      for (double value : {mass, xi0, nu, clover_r, clover_t})
        if (!std::isfinite(value))
          throw std::runtime_error("nonfinite production action parameter");
      if (xi0 <= 0 || nu <= 0 || source_time >= shape[3] ||
          first_time + time_count > shape[3] || !time_count)
        throw std::runtime_error("invalid production action geometry or times");
    }
  };
  inline void write_header(std::ostream &out, const Header &h) {
    h.validate();
    chroma_bridge::write(out, magic.data(), magic.size());
    chroma_bridge::write(out, h.shape.data(), 4);
    for (double value : {h.mass, h.xi0, h.nu, h.clover_r, h.clover_t})
      chroma_bridge::write(out, &value, 1);
    for (uint32_t value : {h.source_time, h.first_time, h.time_count})
      chroma_bridge::write(out, &value, 1);
  }
  inline Header read_header(std::istream &in) {
    Header h;
    std::array<char, 8> tag;
    chroma_bridge::read(in, tag.data(), 8);
    bool aligned =
        tag == std::array<char, 8>{'C', 'H', 'R', 'M', 'G', '0', '0', '3'};
    if (tag != magic && !aligned)
      throw std::runtime_error("unsupported production bridge");
    chroma_bridge::read(in, h.shape.data(), 4);
    for (auto *p : {&h.mass, &h.xi0, &h.nu, &h.clover_r, &h.clover_t})
      chroma_bridge::read(in, p, 1);
    for (auto *p : {&h.source_time, &h.first_time, &h.time_count})
      chroma_bridge::read(in, p, 1);
    if (aligned) {
      uint32_t padding;
      chroma_bridge::read(in, &padding, 1);
      if (padding)
        throw std::runtime_error("bridge padding");
    }
    h.validate();
    return h;
  }
} // namespace production
