#pragma once
#include "wilson_common.hpp"

struct EvenOdd {
  const System &full;
  System eo, oe;
  std::vector<uint32_t> even, odd, index, parity;
  double diagonal;
  EvenOdd(const WilsonSystem &w) : full(w.a), diagonal(w.a.mass + 4) {
    if (w.extent % 2)
      throw std::runtime_error(
          "even/odd requires even periodic lattice extents");
    const auto &p = full.p;
    index.resize(p.rows);
    parity.resize(p.rows);
    for (uint32_t i = 0; i < p.rows; ++i) {
      uint32_t q = i, par = 0;
      for (int mu = 0; mu < 4; ++mu) {
        par += q % w.extent;
        q /= w.extent;
      }
      parity[i] = par % 2;
      auto &sites = parity[i] ? odd : even;
      index[i] = uint32_t(sites.size());
      sites.push_back(i);
    }
    if (even.size() != odd.size())
      throw std::runtime_error("unbalanced checkerboard");
    auto extract = [&](System &dest, const std::vector<uint32_t> &sites) {
      dest.mass = full.mass;
      auto &v = dest.p;
      v.rows = uint32_t(sites.size());
      v.rhs = p.rhs;
      v.colors = 3;
      v.spins = 4;
      v.row.push_back(0);
      for (auto i : sites) {
        if (p.row[i + 1] - p.row[i] != 17 || p.col[p.row[i]] != i)
          throw std::runtime_error("unsupported Wilson diagonal");
        for (uint32_t e = p.row[i] + 1; e < p.row[i + 1]; ++e) {
          if (parity[p.col[e]] == parity[i])
            throw std::runtime_error("hopping edge fails checkerboard parity");
          v.col.push_back(index[p.col[e]]);
          for (int j = 0; j < 9; ++j) {
            dest.ad.push_back(full.ad[e * 9 + j]);
            dest.af.push_back(full.af[e * 9 + j]);
            v.a.push_back(p.a[e * 9 + j]);
          }
          for (int j = 0; j < 4; ++j) {
            v.perm.push_back(p.perm[e * 4 + j]);
            v.phase.push_back(p.phase[e * 4 + j]);
          }
        }
        v.row.push_back(uint32_t(v.col.size()));
      }
      v.x.resize(v.outputs());
    };
    extract(eo, even);
    extract(oe, odd);
  }
};

template <class T> struct CPUSchur {
  const EvenOdd &e;
  int workers;
  std::vector<Cx<T>> temp;
  CPUSchur(const EvenOdd &e, int workers)
      : e(e), workers(workers), temp(e.eo.p.outputs()) {}
  void apply(const Cx<T> *x, Cx<T> *y) {
    apply_cpu(e.oe, x, temp.data(), workers);
    apply_cpu(e.eo, temp.data(), y, workers);
    T inv = T(1 / (e.diagonal * e.diagonal));
    uint32_t rhs = e.eo.p.rhs;
    chunks(e.eo.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j)
        y[j] = {x[j].re - inv * y[j].re, x[j].im - inv * y[j].im};
    });
  }
};

struct GPUSchur {
  GPUProblem eo, oe;
  id<MTLBuffer> temp;
  id<MTLComputePipelineState> combine;
  struct Params {
    uint32_t count;
    float inv_d2;
  } params;
  GPUSchur(MetalRun &m, const EvenOdd &e) : eo(m, e.eo.p), oe(m, e.oe.p) {
    temp = [m.device newBufferWithLength:e.eo.p.outputs() * 8
                                 options:MTLResourceStorageModeShared];
    if (!temp)
      throw std::runtime_error("Schur buffer allocation failed");
    NSError *error = nil;
    auto fn = [m.library newFunctionWithName:@"schur_combine"];
    combine = [m.device newComputePipelineStateWithFunction:fn error:&error];
    if (!combine)
      throw std::runtime_error([[error localizedDescription] UTF8String]);
    params = {uint32_t(e.eo.p.outputs()), float(1 / (e.diagonal * e.diagonal))};
  }
  void encode(id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
              id<MTLBuffer> y) {
    auto hop = [&](GPUProblem &op, id<MTLBuffer> in, id<MTLBuffer> out) {
      [enc setComputePipelineState:op.pipeline];
      for (size_t j = 0; j < op.buffers.size(); ++j)
        [enc setBuffer:op.buffers[j] offset:0 atIndex:j];
      [enc setBuffer:in offset:0 atIndex:3];
      [enc setBuffer:out offset:0 atIndex:4];
      [enc setBytes:&op.rows length:4 atIndex:7];
      [enc dispatchThreads:MTLSizeMake(op.outputs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    };
    hop(oe, x, temp);
    hop(eo, temp, y);
    [enc setComputePipelineState:combine];
    [enc setBuffer:x offset:0 atIndex:0];
    [enc setBuffer:y offset:0 atIndex:1];
    [enc setBytes:&params length:sizeof(params) atIndex:2];
    [enc dispatchThreads:MTLSizeMake(params.count, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(128, 1, 1)];
    [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
  }
};

// Adapter to the existing full-volume FP64 refinement loop. Reduction of the
// RHS and reconstruction are both FP64 CPU work and are included in solve
// timing.
template <class Inner> struct EOCorrection {
  const EvenOdd &e;
  Inner &inner;
  int workers;
  std::vector<D> reduced, odd_rhs, aux, xe, delta;
  double prepare_ms = 0, reconstruct_ms = 0;
  EOCorrection(const EvenOdd &e, Inner &inner, int workers)
      : e(e), inner(inner), workers(workers), reduced(e.eo.p.outputs()),
        odd_rhs(reduced.size()), aux(reduced.size()), xe(reduced.size()),
        delta(e.full.p.outputs()) {}
  void reset_metrics() { prepare_ms = reconstruct_ms = 0; }
  InnerInfo run(const std::vector<D> &residual, double tol, int limit) {
    auto start = Clock::now();
    uint32_t width = e.full.p.rhs * 12;
    chunks(e.eo.p.rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t i = first; i < last; ++i)
        for (uint32_t k = 0; k < width; ++k) {
          reduced[i * width + k] = residual[size_t(e.even[i]) * width + k];
          odd_rhs[i * width + k] = residual[size_t(e.odd[i]) * width + k];
        }
    });
    apply_cpu(e.eo, odd_rhs.data(), aux.data(), workers);
    double inv = 1 / e.diagonal;
    chunks(e.eo.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      uint32_t rhs = e.full.p.rhs;
      for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j)
        reduced[j] = {(reduced[j].re - inv * aux[j].re) * inv,
                      (reduced[j].im - inv * aux[j].im) * inv};
    });
    prepare_ms += seconds(start) * 1000;
    InnerInfo info = inner.run(reduced, tol, limit);
    start = Clock::now();
    auto *de = inner.correction();
    chunks(e.eo.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      uint32_t rhs = e.full.p.rhs;
      for (size_t j = size_t(first) * rhs; j < size_t(last) * rhs; ++j)
        xe[j] = {double(de[j].re), double(de[j].im)};
    });
    apply_cpu(e.oe, xe.data(), aux.data(), workers);
    chunks(e.eo.p.rows, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t i = first; i < last; ++i)
        for (uint32_t k = 0; k < width; ++k) {
          size_t j = i * width + k;
          delta[size_t(e.even[i]) * width + k] = xe[j];
          delta[size_t(e.odd[i]) * width + k] = {
              (odd_rhs[j].re - aux[j].re) * inv,
              (odd_rhs[j].im - aux[j].im) * inv};
        }
    });
    reconstruct_ms += seconds(start) * 1000;
    return info;
  }
  const D *correction() const { return delta.data(); }
};
