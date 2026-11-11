// Tuning-only vector operations. The Wilson/clover operator is unchanged.
struct TuneVectorParams {
  uint n, mode;
  float alpha;
};
kernel void tune_scale_add(device const float2 *x [[buffer(0)]],
                           device float2 *y [[buffer(1)]],
                           constant TuneVectorParams &p [[buffer(2)]],
                           uint i [[thread_position_in_grid]]) {
  if (i >= p.n)
    return;
  if (p.mode == 2)
    y[i] = 0;
  else if (p.mode == 1)
    y[i] = fma(float2(p.alpha), x[i], y[i]);
  else
    y[i] = p.alpha * x[i];
}
