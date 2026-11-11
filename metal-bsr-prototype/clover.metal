struct CloverParams {
  uint rows;
  float beta;
};
kernel void clover_local(device const float2 *a [[buffer(0)]],
                         device const float2 *x [[buffer(1)]],
                         device float2 *y [[buffer(2)]],
                         constant CloverParams &p [[buffer(3)]],
                         uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * 12 * RHS)
    return;
  uint r = gid % RHS, s = (gid / RHS) % 4, c = (gid / (RHS * 4)) % 3,
       site = gid / (12 * RHS), ch = s / 2, q = (s % 2) * 3 + c;
  float2 acc = 0;
  for (uint k = 0; k < 6; ++k) {
    float2 z = a[(site * 2 + ch) * 36 + q * 6 + k],
           v = x[((site * 3 + k % 3) * 4 + ch * 2 + k / 3) * RHS + r];
    acc.x = fma(z.x, v.x, acc.x);
    acc.x = fma(-z.y, v.y, acc.x);
    acc.y = fma(z.x, v.y, acc.y);
    acc.y = fma(z.y, v.x, acc.y);
  }
  if (p.beta != 0)
    acc = fma(float2(p.beta), y[gid], acc);
  y[gid] = acc;
}
kernel void clover_subtract(device const float2 *a [[buffer(0)]],
                            device const float2 *b [[buffer(1)]],
                            device float2 *y [[buffer(2)]],
                            constant CloverParams &p [[buffer(3)]],
                            uint gid [[thread_position_in_grid]]) {
  if (gid < p.rows * 12 * RHS)
    y[gid] = a[gid] - b[gid];
}
kernel void clover_gather(device const float2 *full [[buffer(0)]],
                          device const uint *even [[buffer(1)]],
                          device const uint *odd [[buffer(2)]],
                          device float2 *re [[buffer(3)]],
                          device float2 *ro [[buffer(4)]],
                          constant CloverParams &p [[buffer(5)]],
                          uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * 12 * RHS)
    return;
  uint site = gid / (12 * RHS), q = gid % (12 * RHS);
  re[gid] = full[even[site] * 12 * RHS + q];
  ro[gid] = full[odd[site] * 12 * RHS + q];
}
kernel void clover_scatter(device const float2 *de [[buffer(0)]],
                           device const float2 *dodd [[buffer(1)]],
                           device float2 *full [[buffer(2)]],
                           device const uint *even [[buffer(3)]],
                           device const uint *odd [[buffer(4)]],
                           constant CloverParams &p [[buffer(5)]],
                           uint gid [[thread_position_in_grid]]) {
  if (gid >= p.rows * 12 * RHS)
    return;
  uint site = gid / (12 * RHS), q = gid % (12 * RHS);
  full[even[site] * 12 * RHS + q] = de[gid];
  full[odd[site] * 12 * RHS + q] = dodd[gid];
}
