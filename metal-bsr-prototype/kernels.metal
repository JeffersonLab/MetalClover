#include <metal_stdlib>
using namespace metal;

// Specialize shapes, not numerical data. Both variants have 12 output components.
constant uint RHS [[function_constant(0)]];
constant uint COLORS [[function_constant(1)]];
constant uint SPINS [[function_constant(2)]];

struct Params { uint rows; };

inline float2 cmul(float2 a, float2 b) {
    return float2(fma(a.x, b.x, -a.y * b.y), fma(a.x, b.y, a.y * b.x));
}

// Generic CSR of dense blocks, or compressed A_color (x) spin permutation.
// Every thread owns one complex output: no atomics or inter-thread reductions.
kernel void bsr_multiply(device const uint* row_ptr [[buffer(0)]],
                         device const uint* col_idx [[buffer(1)]],
                         device const float2* a [[buffer(2)]],
                         device const float2* x [[buffer(3)]],
                         device float2* y [[buffer(4)]],
                         device const uint* perm [[buffer(5)]],
                         device const float2* phase [[buffer(6)]],
                         constant Params& p [[buffer(7)]],
                         uint gid [[thread_position_in_grid]]) {
    if (gid >= p.rows * COLORS * SPINS * RHS) return;
    uint r = gid % RHS;
    uint s = (gid / RHS) % SPINS;
    uint c = (gid / (RHS * SPINS)) % COLORS;
    uint row = gid / (RHS * SPINS * COLORS);
    float2 acc = float2(0);
    for (uint e = row_ptr[row]; e < row_ptr[row + 1]; ++e) {
        uint si = SPINS == 1 ? 0 : perm[e * SPINS + s];
        float2 ph = SPINS == 1 ? float2(1, 0) : phase[e * SPINS + s];
        for (uint k = 0; k < COLORS; ++k) {
            float2 av = a[(e * COLORS + c) * COLORS + k];
            float2 xv = x[((col_idx[e] * COLORS + k) * SPINS + si) * RHS + r];
            if (SPINS != 1) xv = cmul(ph, xv);
            acc.x = fma(av.x, xv.x, acc.x);
            acc.x = fma(-av.y, xv.y, acc.x);
            acc.y = fma(av.x, xv.y, acc.y);
            acc.y = fma(av.y, xv.x, acc.y);
        }
    }
    y[gid] = acc;
}
