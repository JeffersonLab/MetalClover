#pragma once
#include "clover_common.hpp"

uint64_t rhs_hash(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
template <class Clover>
void set_rhs(Clover &c, uint32_t offset, int workers, bool zero = false) {
  auto &a = c.w.a;
  uint32_t rhs = a.p.rhs;
  chunks(a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
    for (size_t q = first; q < last; ++q)
      for (uint32_t r = 0; r < rhs; ++r) {
        uint64_t seed = rhs_hash(q) ^ rhs_hash(uint64_t(offset + r) + 771337);
        auto uniform = [](uint64_t v) {
          return 2 * double(v >> 11) * 0x1p-53 - 1;
        };
        a.truth[q * rhs + r] =
            zero && r == rhs - 1
                ? D{}
                : D{uniform(rhs_hash(seed)), uniform(rhs_hash(seed + 1))};
      }
  });
  apply_cpu(a, a.truth.data(), a.b.data(), workers);
  a.bnorm.assign(rhs, 0);
  std::vector<double> partial(size_t(workers) * rhs);
  chunks(a.p.rows * 12, workers, [&](int w, uint32_t first, uint32_t last) {
    for (size_t q = first; q < last; ++q)
      for (uint32_t r = 0; r < rhs; ++r) {
        auto z = a.b[q * rhs + r];
        partial[size_t(w) * rhs + r] += z.re * z.re + z.im * z.im;
      }
  });
  for (uint32_t r = 0; r < rhs; ++r) {
    for (int w = 0; w < workers; ++w)
      a.bnorm[r] += partial[size_t(w) * rhs + r];
    a.bnorm[r] = std::sqrt(a.bnorm[r]);
  }
}
