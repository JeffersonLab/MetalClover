struct SchurParams { uint count; float inv_d2; };
kernel void schur_combine(device const float2* x [[buffer(0)]],
                          device float2* y [[buffer(1)]],
                          constant SchurParams& p [[buffer(2)]],
                          uint gid [[thread_position_in_grid]]) {
    if(gid<p.count)y[gid]=fma(float2(-p.inv_d2),y[gid],x[gid]);
}
