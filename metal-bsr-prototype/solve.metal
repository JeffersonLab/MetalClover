// Appended to kernels.metal at runtime: RHS is the same function constant.
struct CGParams {
  uint components, parts, max_iterations;
  float tolerance2;
};
// scalars: rho, initial_rho, alpha, beta, active, iterations, breakdown.
kernel void cg_partial(device const float2 *v [[buffer(0)]],
                       device const float2 *w [[buffer(1)]],
                       device float *partial [[buffer(2)]],
                       constant CGParams &p [[buffer(3)]],
                       uint gid [[thread_position_in_grid]]) {
  if (gid >= RHS * p.parts)
    return;
  uint r = gid % RHS, part = gid / RHS;
  float sum = 0;
  for (uint j = part * 256; j < min((part + 1) * 256, p.components); ++j) {
    float2 a = v[j * RHS + r], b = w[j * RHS + r];
    sum = fma(a.x, b.x, sum);
    sum = fma(a.y, b.y, sum);
  }
  partial[gid] = sum;
}
kernel void cg_init(device const float *partial [[buffer(0)]],
                    device float *scalars [[buffer(1)]],
                    constant CGParams &p [[buffer(2)]],
                    uint r [[thread_position_in_grid]]) {
  if (r >= RHS)
    return;
  float rho = 0;
  for (uint k = 0; k < p.parts; ++k)
    rho += partial[k * RHS + r];
  for (uint j = 0; j < 7; ++j)
    scalars[j * RHS + r] = 0;
  scalars[r] = rho;
  scalars[RHS + r] = rho;
  scalars[4 * RHS + r] = rho > 0 ? 1 : 0;
}
kernel void cg_alpha(device const float *partial [[buffer(0)]],
                     device float *scalars [[buffer(1)]],
                     constant CGParams &p [[buffer(2)]],
                     uint r [[thread_position_in_grid]]) {
  if (r >= RHS)
    return;
  float pap = 0;
  for (uint k = 0; k < p.parts; ++k)
    pap += partial[k * RHS + r];
  float alpha = 0;
  if (scalars[4 * RHS + r] > 0) {
    if (!(pap > 0) || !isfinite(pap)) {
      scalars[6 * RHS + r] = 1;
      scalars[4 * RHS + r] = 0;
    } else
      alpha = scalars[r] / pap;
  }
  scalars[2 * RHS + r] = alpha;
}
kernel void cg_update(device float2 *delta [[buffer(0)]],
                      device float2 *residual [[buffer(1)]],
                      device const float2 *direction [[buffer(2)]],
                      device const float2 *ap [[buffer(3)]],
                      device const float *scalars [[buffer(4)]],
                      constant CGParams &p [[buffer(5)]],
                      uint gid [[thread_position_in_grid]]) {
  if (gid >= p.components * RHS)
    return;
  float alpha = scalars[2 * RHS + gid % RHS];
  delta[gid] = fma(float2(alpha), direction[gid], delta[gid]);
  residual[gid] = fma(float2(-alpha), ap[gid], residual[gid]);
}
kernel void cg_beta(device const float *partial [[buffer(0)]],
                    device float *scalars [[buffer(1)]],
                    constant CGParams &p [[buffer(2)]],
                    uint r [[thread_position_in_grid]]) {
  if (r >= RHS)
    return;
  float next = 0;
  for (uint k = 0; k < p.parts; ++k)
    next += partial[k * RHS + r];
  float beta = 0;
  if (scalars[4 * RHS + r] > 0) {
    scalars[5 * RHS + r] += 1;
    if (!isfinite(next)) {
      scalars[6 * RHS + r] = 1;
      scalars[4 * RHS + r] = 0;
    } else if (next <= p.tolerance2 * scalars[RHS + r] ||
               scalars[5 * RHS + r] >= p.max_iterations)
      scalars[4 * RHS + r] = 0;
    else
      beta = next / scalars[r];
  }
  scalars[r] = next;
  scalars[3 * RHS + r] = beta;
}
kernel void cg_direction(device const float2 *residual [[buffer(0)]],
                         device float2 *direction [[buffer(1)]],
                         device const float *scalars [[buffer(2)]],
                         constant CGParams &p [[buffer(3)]],
                         uint gid [[thread_position_in_grid]]) {
  if (gid >= p.components * RHS)
    return;
  uint r = gid % RHS;
  direction[gid] = scalars[4 * RHS + r] > 0 ? fma(float2(scalars[3 * RHS + r]),
                                                  direction[gid], residual[gid])
                                            : float2(0);
}
