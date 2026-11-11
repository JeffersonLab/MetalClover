struct CompactParams {
  uint rows, lattice;
};
kernel void compact_hop(device const float2 *links [[buffer(0)]],
                        device const uint *neighbors [[buffer(1)]],
                        device const uint *map [[buffer(2)]],
                        device const uint *perm [[buffer(3)]],
                        device const float2 *phase [[buffer(4)]],
                        device const float2 *x [[buffer(5)]],
                        device float2 *y [[buffer(6)]],
                        constant CompactParams &p [[buffer(7)]],
                        uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * 12 * RHS)
    return;
  uint r = gid % RHS, q = gid / RHS, s = q % 4, c = q / 4 % 3,
       site = map[q / 12];
  float2 sum = 0;
  for (uint d = 0; d < 8; ++d) {
    uint mu = d / 2, n = neighbors[site * 8 + d], si = perm[mu * 4 + s];
    bool f = d % 2 == 0;
    float2 ph = (f ? -1.f : 1.f) * phase[mu * 4 + s];
    uint time = site / (p.lattice * p.lattice * p.lattice);
    float scale = mu == 3 && ((f && time == p.lattice - 1) || (!f && time == 0))
                      ? .5f
                      : -.5f;
    for (uint k = 0; k < 3; ++k) {
      float2 a =
          links[((f ? site : n) * 4 + mu) * 9 + (f ? c * 3 + k : k * 3 + c)];
      if (!f)
        a.y = -a.y;
      a *= scale;
      float2 v = x[((n / 2) * 12 + k * 4 + s) * RHS + r] +
                 cmul(ph, x[((n / 2) * 12 + k * 4 + si) * RHS + r]);
      sum += cmul(a, v);
    }
  }
  y[gid] = sum;
}
struct MG3Params {
  uint rows, inbs, nvec, coarse, mode, spin;
  float omega;
};
kernel void mg3_diagonal(device const float2 *a [[buffer(0)]],
                         device const float2 *x [[buffer(1)]],
                         device float2 *y [[buffer(2)]],
                         constant MG3Params &p [[buffer(3)]],
                         uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * p.inbs * RHS)
    return;
  uint r = gid % RHS, q = gid / RHS, site = q / p.inbs, i = q % p.inbs;
  float2 sum = 0;
  if (p.spin) {
    uint ch = i % 4 / 2, row = (i % 2) * 3 + i / 4;
    for (uint k = 0; k < 6; ++k)
      sum += cmul(a[(site * 2 + ch) * 36 + row * 6 + k],
                  x[(site * 12 + (k % 3) * 4 + ch * 2 + k / 3) * RHS + r]);
  } else
    for (uint k = 0; k < p.inbs; ++k)
      sum += cmul(a[(site * p.inbs + i) * p.inbs + k],
                  x[(site * p.inbs + k) * RHS + r]);
  y[gid] = sum;
}
kernel void mg3_restrict(device const float2 *x [[buffer(0)]],
                         device const float2 *basis [[buffer(1)]],
                         device const uint *offsets [[buffer(2)]],
                         device const uint *sites [[buffer(3)]],
                         device float2 *y [[buffer(4)]],
                         constant MG3Params &p [[buffer(5)]],
                         uint gid [[thread_position_in_grid]]) {
  if (gid >= p.coarse * RHS)
    return;
  uint r = gid % RHS, row = gid / RHS, ag = row / (2 * p.nvec),
       ch = row % (2 * p.nvec) / p.nvec, v = row % p.nvec;
  float2 sum = 0;
  for (uint j = offsets[ag]; j < offsets[ag + 1]; ++j)
    for (uint q = 0; q < p.inbs; ++q) {
      uint c = p.spin ? q % 4 / 2 : q / (p.inbs / 2);
      if (c == ch) {
        uint i = sites[j] * p.inbs + q;
        float2 w = basis[i * p.nvec + v], z = x[i * RHS + r];
        w.y = -w.y;
        sum += cmul(w, z);
      }
    }
  y[gid] = sum;
}
kernel void mg3_prolong(device const float2 *x [[buffer(0)]],
                        device const float2 *basis [[buffer(1)]],
                        device const uint *aggregate [[buffer(2)]],
                        device float2 *y [[buffer(3)]],
                        constant MG3Params &p [[buffer(4)]],
                        uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * p.inbs * RHS)
    return;
  uint q = gid / RHS, r = gid % RHS,
       ch = p.spin ? q % 4 / 2 : (q % p.inbs) / (p.inbs / 2),
       base = aggregate[q / p.inbs] * 2 * p.nvec + ch * p.nvec;
  float2 sum = 0;
  for (uint v = 0; v < p.nvec; ++v)
    sum += cmul(basis[q * p.nvec + v], x[(base + v) * RHS + r]);
  y[gid] += sum;
}
kernel void mg3_difference(device const float2 *a [[buffer(0)]],
                           device const float2 *b [[buffer(1)]],
                           device float2 *y [[buffer(2)]],
                           constant MG3Params &p [[buffer(3)]],
                           uint gid [[thread_position_in_grid]]) {
  if (gid < p.rows * p.inbs * RHS)
    y[gid] = a[gid] - b[gid];
}
kernel void mg3_add(device const float2 *x [[buffer(0)]],
                    device float2 *y [[buffer(1)]],
                    constant MG3Params &p [[buffer(2)]],
                    uint gid [[thread_position_in_grid]]) {
  if (gid < p.rows * p.inbs * RHS)
    y[gid] += x[gid];
}
kernel void mg3_copy(device const float2 *x [[buffer(0)]],
                     device float2 *y [[buffer(1)]],
                     constant MG3Params &p [[buffer(2)]],
                     uint gid [[thread_position_in_grid]]) {
  y[gid] = x[gid];
}

struct CompactRectParams {
  uint rows, time_stride, time_extent, antiperiodic;
};
kernel void compact_hop_rect(device const float2 *links [[buffer(0)]],
                             device const uint *neighbors [[buffer(1)]],
                             device const uint *map [[buffer(2)]],
                             device const uint *perm [[buffer(3)]],
                             device const float2 *phase [[buffer(4)]],
                             device const float2 *x [[buffer(5)]],
                             device float2 *y [[buffer(6)]],
                             constant CompactRectParams &p [[buffer(7)]],
                             uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * 12 * RHS)
    return;
  uint r = gid % RHS, q = gid / RHS, s = q % 4, c = q / 4 % 3,
       site = map[q / 12];
  float2 sum = 0;
  for (uint d = 0; d < 8; ++d) {
    uint mu = d / 2, n = neighbors[site * 8 + d], si = perm[mu * 4 + s];
    bool f = d % 2 == 0;
    float2 ph = (f ? -1.f : 1.f) * phase[mu * 4 + s];
    uint time = site / p.time_stride;
    float scale =
        p.antiperiodic && mu == 3 &&
                ((f && time == p.time_extent - 1) || (!f && time == 0))
            ? .5f
            : -.5f;
    for (uint k = 0; k < 3; ++k) {
      float2 a =
          links[((f ? site : n) * 4 + mu) * 9 + (f ? c * 3 + k : k * 3 + c)];
      if (!f)
        a.y = -a.y;
      a *= scale;
      float2 v = x[((n / 2) * 12 + k * 4 + s) * RHS + r] +
                 cmul(ph, x[((n / 2) * 12 + k * 4 + si) * RHS + r]);
      sum += cmul(a, v);
    }
  }
  y[gid] = sum;
}
