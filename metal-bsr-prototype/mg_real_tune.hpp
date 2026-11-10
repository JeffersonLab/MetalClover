#pragma once
#include "three_level.hpp"

struct TuneVector {
  MetalRun &m;
  id<MTLComputePipelineState> pipeline;
  TuneVector(MetalRun &m) : m(m) {
    NSError *e = nil;
    auto f = [m.library newFunctionWithName:@"tune_scale_add"];
    pipeline = [m.device newComputePipelineStateWithFunction:f error:&e];
    if (!pipeline)
      throw std::runtime_error([[e localizedDescription] UTF8String]);
  }
  void encode(id<MTLComputeCommandEncoder> e, id<MTLBuffer> x, id<MTLBuffer> y,
              uint32_t n, float alpha, uint32_t mode) {
    struct P {
      uint32_t n, mode;
      float alpha;
    } p{n, mode, alpha};
    [e setComputePipelineState:pipeline];
    [e setBuffer:x offset:0 atIndex:0];
    [e setBuffer:y offset:0 atIndex:1];
    [e setBytes:&p length:sizeof(p) atIndex:2];
    [e dispatchThreads:MTLSizeMake(n, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [e memoryBarrierWithScope:MTLBarrierScopeBuffers];
  }
};

// Defaults reproduce the original one-pre/one-post cycle, including its
// smoother weight 1. Changing cycles does not change the fine operator.
struct TunedGPUCycle : GPUCycle {
  TuneVector vectors;
  int pre_steps = 1, post_steps = 1;
  float omega = 1;
  TunedGPUCycle(const Transfer &h, MetalRun &m, uint32_t rhs,
                const std::vector<Cx<float>> &inverse, uint32_t format = 0)
      : GPUCycle(h, m, rhs, inverse, format), vectors(m) {}
  void diagonal(id<MTLComputeCommandEncoder> e, id<MTLBuffer> r,
                id<MTLBuffer> z) {
    size_t n = size_t(h.rows) * h.inbs * rhs;
    dispatch(e, 0, {inv, r, z}, n);
    if (omega != 1)
      vectors.encode(e, z, z, uint32_t(n), omega, 0);
  }
  void residual_for(id<MTLComputeCommandEncoder> e, id<MTLBuffer> r,
                    id<MTLBuffer> z) {
    product(e, z, work);
    ++products;
    dispatch(e, 3, {r, work, residual}, size_t(h.rows) * h.inbs * rhs);
  }
  void before(id<MTLComputeCommandEncoder> e, id<MTLBuffer> r,
              id<MTLBuffer> z) {
    ++calls;
    size_t n = size_t(h.rows) * h.inbs * rhs;
    if (pre_steps) {
      diagonal(e, r, z);
      for (int i = 1; i < pre_steps; ++i) {
        residual_for(e, r, z);
        diagonal(e, residual, smooth);
        dispatch(e, 4, {smooth, z}, n);
      }
    } else
      vectors.encode(e, r, z, uint32_t(n), 0, 2);
    if (enabled) {
      if (pre_steps)
        residual_for(e, r, z);
      else
        dispatch(e, 5, {r, residual}, n);
      dispatch(e, 1, {residual, basis, offsets, sites, b},
               size_t(h.coarse) * rhs);
    }
  }
  void after(id<MTLComputeCommandEncoder> e, id<MTLBuffer> r, id<MTLBuffer> z) {
    size_t n = size_t(h.rows) * h.inbs * rhs;
    if (enabled)
      dispatch(e, 2, {x, basis, aggregate, z}, n);
    for (int i = 0; i < post_steps; ++i) {
      if (enabled || pre_steps || i)
        residual_for(e, r, z);
      else
        dispatch(e, 5, {r, residual}, n);
      diagonal(e, residual, smooth);
      dispatch(e, 4, {smooth, z}, n);
    }
  }
  void encode(id<MTLComputeCommandEncoder> e, id<MTLBuffer> r,
              id<MTLBuffer> z) {
    before(e, r, z);
    if (enabled)
      coarse(e, b, x);
    after(e, r, z);
  }
  void apply_device(id<MTLBuffer> r, id<MTLBuffer> z) {
    @autoreleasepool {
      auto cmd = [m.queue commandBuffer];
      auto e = [cmd computeCommandEncoder];
      before(e, r, z);
      finish(cmd, e);
      if (enabled)
        blocking_coarse(b, x);
      cmd = [m.queue commandBuffer];
      e = [cmd computeCommandEncoder];
      after(e, r, z);
      finish(cmd, e);
    }
  }
};

// Filter the same deterministic starting vectors as Transfer::build, on the
// GPU in FP32. Local QR and sparse Galerkin assembly remain CPU FP64.
// The callback writes x-filter(x), then the existing build update with omega=1
// installs filter(x). Its one-vector CPU tile bounds setup memory.
template <class Schur>
void build_gpu_filtered_basis(Transfer &h, Schur &op, const CompactClover &c,
                              int sweeps, double omega, int workers) {
  TuneVector vectors(op.m);
  size_t n = size_t(h.rows) * h.inbs;
  auto alloc = [&](const void *p, size_t bytes) {
    auto b = [op.m.device newBufferWithLength:bytes
                                      options:MTLResourceStorageModeShared];
    if (!b)
      throw std::runtime_error("GPU basis filter allocation");
    if (p)
      std::memcpy([b contents], p, bytes);
    return b;
  };
  auto inverse = c.factorized
                     ? op.Aodd
                     : alloc(c.inverseOdd32.data(), c.inverseOdd32.size() * 8),
       x = alloc(nullptr, n * 8), y = alloc(nullptr, n * 8);
  h.build(
      [&](const D *input, D *difference) {
        auto *values = static_cast<Cx<float> *>([x contents]);
        chunks(uint32_t(n), workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = first; i < last; ++i)
            values[i] = {float(input[i].re), float(input[i].im)};
        });
        for (int first = 0; first < sweeps; first += 8) {
          @autoreleasepool {
            auto cmd = [op.m.queue commandBuffer];
            auto e = [cmd computeCommandEncoder];
            for (int k = first; k < std::min(sweeps, first + 8); ++k) {
              op.encode(e, x, y);
              op.block(e, inverse, y, op.inverted);
              vectors.encode(e, op.inverted, x, uint32_t(n), float(-omega), 1);
            }
            [e endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
            if (cmd.status == MTLCommandBufferStatusError)
              throw std::runtime_error(
                  [[cmd.error localizedDescription] UTF8String]);
          }
        }
        chunks(uint32_t(n), workers, [&](int, uint32_t first, uint32_t last) {
          for (size_t i = first; i < last; ++i)
            difference[i] = {input[i].re - double(values[i].re),
                             input[i].im - double(values[i].im)};
        });
      },
      [&](const D *x, D *y) { std::copy_n(x, n, y); }, 1, 1, 1);
}
