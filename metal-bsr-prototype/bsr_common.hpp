#pragma once
#include "operator_work.hpp"
#include <Accelerate/Accelerate.h>
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <algorithm>
#include <arm_neon.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <dispatch/dispatch.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MetalCloverEnv {

  using Clock = std::chrono::steady_clock;
  struct Z {
    float re, im;
  };
  static_assert(sizeof(Z) == 8 && sizeof(std::complex<float>) == 8);
  double seconds(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
  }
  Z mul(Z a, Z b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
  }
  Z add(Z a, Z b) { return {a.re + b.re, a.im + b.im}; }
  const __LAPACK_float_complex *blas_ptr(const Z *p) {
    return reinterpret_cast<const __LAPACK_float_complex *>(p);
  }
  __LAPACK_float_complex *blas_ptr(Z *p) {
    return reinterpret_cast<__LAPACK_float_complex *>(p);
  }

  struct Problem {
    uint32_t rows, rhs, colors, spins;
    std::vector<uint32_t> row, col, perm;
    std::vector<Z> a, phase, x;
    std::string mode() const { return spins == 1 ? "dense12" : "color3_spin4"; }
    size_t outputs() const { return size_t(rows) * colors * spins * rhs; }
    double flops() const {
      // Useful work: 8 real FLOPs per color complex multiply-add. Exclude
      // phase/packing.
      return double(col.size()) * colors * spins * rhs * 8 * colors;
    }
    size_t bytes() const {
      return (row.size() + col.size() + perm.size()) * 4 +
             (a.size() + phase.size() + x.size() + outputs()) * 8;
    }
  };

  Problem make_problem(uint32_t extent, uint32_t rhs, bool kron,
                       bool irregular = false) {
    Problem p;
    p.rows = extent * extent * extent * extent;
    p.rhs = rhs;
    p.colors = kron ? 3 : 12;
    p.spins = kron ? 4 : 1;
    std::mt19937 gen(20261103);
    std::uniform_real_distribution<float> dist(-0.25f, 0.25f);
    const Z phases[] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    p.row.push_back(0);
    for (uint32_t i = 0; i < p.rows; ++i) {
      for (uint32_t d = 0; d < 9; ++d) {
        if (irregular && (i % 7 == 0 || (i + d) % 5 == 0))
          continue;
        uint32_t j = i;
        if (d > 0) {
          uint32_t stride = 1;
          for (uint32_t mu = 0; mu < (d - 1) / 2; ++mu)
            stride *= extent;
          uint32_t coord = (i / stride) % extent;
          int delta = d % 2 ? (coord + 1 == extent ? -int(extent - 1) : 1)
                            : (coord == 0 ? int(extent - 1) : -1);
          j = uint32_t(int64_t(i) + int64_t(delta) * stride);
        }
        p.col.push_back(j);
        for (uint32_t c = 0; c < p.colors * p.colors; ++c)
          p.a.push_back({dist(gen), dist(gen)});
        // Per-edge permutation and unit phases, not an actual gauge
        // configuration.
        for (uint32_t s = 0; s < p.spins; ++s) {
          p.perm.push_back((s + d) % p.spins);
          p.phase.push_back(phases[(s + d) % 4]);
        }
      }
      p.row.push_back(uint32_t(p.col.size()));
    }
    p.x.resize(p.outputs());
    for (auto &z : p.x)
      z = {dist(gen), dist(gen)};
    return p;
  }

  // Keep work partitioning explicit: GCD workers are not an OpenMP dependency.
  void parallel_rows(uint32_t rows, int workers,
                     const std::function<void(uint32_t, uint32_t)> &f) {
    if (workers == 1) {
      f(0, rows);
      return;
    }
    dispatch_apply(size_t(workers),
                   dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t w) {
                     f(uint32_t(uint64_t(rows) * w / workers),
                       uint32_t(uint64_t(rows) * (w + 1) / workers));
                   });
  }

  // ARM NEON: vectorize over four RHS values, retaining output accumulators in
  // registers.
  void cpu_neon(const Problem &p, Z *y, int workers) {
    parallel_rows(p.rows, workers, [&](uint32_t begin, uint32_t end) {
      for (uint32_t i = begin; i < end; ++i)
        for (uint32_t c = 0; c < p.colors; ++c)
          for (uint32_t s = 0; s < p.spins; ++s) {
            uint32_t r = 0;
            for (; r + 4 <= p.rhs; r += 4) {
              float32x4_t re = vdupq_n_f32(0), im = vdupq_n_f32(0);
              for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e) {
                uint32_t si = p.spins == 1 ? 0 : p.perm[e * p.spins + s];
                Z ph = p.spins == 1 ? Z{1, 0} : p.phase[e * p.spins + s];
                for (uint32_t k = 0; k < p.colors; ++k) {
                  Z av = p.a[(size_t(e) * p.colors + c) * p.colors + k];
                  auto xv = vld2q_f32(reinterpret_cast<const float *>(
                      &p.x[((size_t(p.col[e]) * p.colors + k) * p.spins + si) *
                               p.rhs +
                           r]));
                  auto xr = xv.val[0], xi = xv.val[1];
                  if (p.spins != 1) {
                    xr = vfmsq_n_f32(vmulq_n_f32(xv.val[0], ph.re), xv.val[1],
                                     ph.im);
                    xi = vfmaq_n_f32(vmulq_n_f32(xv.val[1], ph.re), xv.val[0],
                                     ph.im);
                  }
                  re = vfmaq_n_f32(re, xr, av.re);
                  re = vfmsq_n_f32(re, xi, av.im);
                  im = vfmaq_n_f32(im, xi, av.re);
                  im = vfmaq_n_f32(im, xr, av.im);
                }
              }
              float32x4x2_t out = {{re, im}};
              vst2q_f32(
                  reinterpret_cast<float *>(
                      &y[((size_t(i) * p.colors + c) * p.spins + s) * p.rhs +
                         r]),
                  out);
            }
            for (; r < p.rhs; ++r) {
              Z acc{0, 0};
              for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e) {
                uint32_t si = p.spins == 1 ? 0 : p.perm[e * p.spins + s];
                Z ph = p.spins == 1 ? Z{1, 0} : p.phase[e * p.spins + s];
                for (uint32_t k = 0; k < p.colors; ++k) {
                  Z xv =
                      p.x[((size_t(p.col[e]) * p.colors + k) * p.spins + si) *
                              p.rhs +
                          r];
                  if (p.spins != 1)
                    xv = mul(ph, xv);
                  acc = add(
                      acc,
                      mul(p.a[(size_t(e) * p.colors + c) * p.colors + k], xv));
                }
              }
              y[((size_t(i) * p.colors + c) * p.spins + s) * p.rhs + r] = acc;
            }
          }
    });
  }

  // Mirrors the per-nonzero-block BLAS strategy in SuperbBLAS's non-MKL CPU
  // path. Spin mode gathers/permutates input before a 3 x (4*RHS) complex GEMM.
  // Scratch is preallocated, but its population and output zeroing are timed.
  void cpu_blas(const Problem &p, Z *y, int workers,
                std::vector<std::vector<Z>> &scratch) {
    parallel_rows(p.rows, workers, [&](uint32_t begin, uint32_t end) {
      uint32_t worker = 0;
      // Recover the partition index even when rows < workers.
      for (int w = 0; w < workers; ++w)
        if (uint32_t(uint64_t(p.rows) * w / workers) == begin &&
            uint32_t(uint64_t(p.rows) * (w + 1) / workers) == end) {
          worker = w;
          break;
        }
      auto &aux = scratch[worker];
      Z one{1, 0};
      for (uint32_t i = begin; i < end; ++i) {
        Z *yi = y + size_t(i) * p.colors * p.spins * p.rhs;
        std::fill(yi, yi + p.colors * p.spins * p.rhs, Z{0, 0});
        for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e) {
          const Z *xi =
              p.x.data() + size_t(p.col[e]) * p.colors * p.spins * p.rhs;
          if (p.spins != 1) {
            for (uint32_t k = 0; k < p.colors; ++k)
              for (uint32_t s = 0; s < p.spins; ++s)
                for (uint32_t r = 0; r < p.rhs; ++r)
                  aux[(k * p.spins + s) * p.rhs + r] = mul(
                      p.phase[e * p.spins + s],
                      xi[(k * p.spins + p.perm[e * p.spins + s]) * p.rhs + r]);
            xi = aux.data();
          }
          if (p.spins == 1 && p.rhs == 1) {
            cblas_cgemv(CblasRowMajor, CblasNoTrans, int(p.colors),
                        int(p.colors), blas_ptr(&one),
                        blas_ptr(p.a.data() + size_t(e) * p.colors * p.colors),
                        int(p.colors), blas_ptr(xi), 1, blas_ptr(&one),
                        blas_ptr(yi), 1);
          } else {
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        int(p.colors), int(p.spins * p.rhs), int(p.colors),
                        blas_ptr(&one),
                        blas_ptr(p.a.data() + size_t(e) * p.colors * p.colors),
                        int(p.colors), blas_ptr(xi), int(p.spins * p.rhs),
                        blas_ptr(&one), blas_ptr(yi), int(p.spins * p.rhs));
          }
        }
      }
    });
  }

  // Independent double-precision scalar oracle, with a different accumulation
  // grouping.
  std::complex<double> reference_one(const Problem &p, size_t idx) {
    uint32_t r = idx % p.rhs, s = (idx / p.rhs) % p.spins,
             c = (idx / (p.rhs * p.spins)) % p.colors;
    uint32_t i = uint32_t(idx / (p.rhs * p.spins * p.colors));
    std::complex<double> sum = 0;
    for (uint32_t e = p.row[i]; e < p.row[i + 1]; ++e) {
      uint32_t si = p.spins == 1 ? 0 : p.perm[e * p.spins + s];
      Z ph = p.spins == 1 ? Z{1, 0} : p.phase[e * p.spins + s];
      std::complex<double> inner = 0;
      for (uint32_t k = 0; k < p.colors; ++k) {
        auto a = p.a[(size_t(e) * p.colors + c) * p.colors + k];
        auto x =
            p.x[((size_t(p.col[e]) * p.colors + k) * p.spins + si) * p.rhs + r];
        inner +=
            std::complex<double>(a.re, a.im) * std::complex<double>(x.re, x.im);
      }
      sum += inner * std::complex<double>(ph.re, ph.im);
    }
    return sum;
  }

  struct Check {
    size_t count = 0;
    double max_abs = 0, rel_l2 = 0;
  };
  std::vector<size_t> check_indices(const Problem &p, bool full) {
    std::vector<size_t> v;
    if (full || p.outputs() <= 65536) {
      v.resize(p.outputs());
      for (size_t j = 0; j < v.size(); ++j)
        v[j] = j;
    } else {
      std::mt19937_64 rng(7341);
      for (size_t j = 0; j < 4096; ++j)
        v.push_back(rng() % p.outputs());
      // Explicitly cover ends, RHS/component boundaries and partition
      // boundaries.
      for (size_t i :
           {size_t(0), size_t(1), size_t(p.rows / 2), size_t(p.rows - 1)})
        for (size_t j = 0; j < p.colors * p.spins * p.rhs; ++j)
          v.push_back(i * p.colors * p.spins * p.rhs + j);
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return v;
  }
  Check validate(const Z *y, const std::vector<size_t> &indices,
                 const std::vector<std::complex<double>> &ref) {
    Check c;
    double err2 = 0, ref2 = 0;
    for (size_t j = 0; j < indices.size(); ++j) {
      Z z = y[indices[j]];
      if (!std::isfinite(z.re) || !std::isfinite(z.im))
        throw std::runtime_error("nonfinite result");
      double err = std::abs(std::complex<double>(z.re, z.im) - ref[j]);
      c.max_abs = std::max(c.max_abs, err);
      err2 += err * err;
      ref2 += std::norm(ref[j]);
      ++c.count;
    }
    c.rel_l2 = std::sqrt(err2 / std::max(ref2, 1e-300));
    if (c.rel_l2 > 5e-6 || c.max_abs > 2e-5)
      throw std::runtime_error("correctness check failed");
    return c;
  }

  struct MetalRun {
    OperatorWork work;
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLLibrary> library;
    MetalRun(const std::string &source) {
      device = MTLCreateSystemDefaultDevice();
      if (!device)
        throw std::runtime_error("No Metal device visible");
      queue = [device newCommandQueue];
      NSError *error = nil;
      auto opts = [MTLCompileOptions new];
      // The main executable may be GCC-built and advertise an older SDK.
      // Do not let its default select Metal 1.1 when loading this dylib.
      opts.languageVersion = MTLLanguageVersion3_1;
      if (@available(macOS 15.0, *))
        opts.mathMode = MTLMathModeSafe;
      else {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        opts.fastMathEnabled = NO;
#pragma clang diagnostic pop
      }
      library = [device
          newLibraryWithSource:[NSString stringWithUTF8String:source.c_str()]
                       options:opts
                         error:&error];
      if (!library)
        throw std::runtime_error([[error localizedDescription] UTF8String]);
    }
  };
  struct GPUProblem {
    MetalRun &m;
    id<MTLComputePipelineState> pipeline;
    std::vector<id<MTLBuffer>> buffers;
    uint32_t rows;
    size_t outputs, threads;
    double upload_ms = 0, compile_ms = 0;
    GPUProblem(MetalRun &m, const Problem &p)
        : m(m), rows(p.rows), outputs(p.outputs()) {
      auto t = Clock::now();
      NSError *error = nil;
      auto values = [MTLFunctionConstantValues new];
      uint32_t rhs = p.rhs, c = p.colors, s = p.spins;
      [values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
      [values setConstantValue:&c type:MTLDataTypeUInt atIndex:1];
      [values setConstantValue:&s type:MTLDataTypeUInt atIndex:2];
      auto fn = [m.library newFunctionWithName:@"bsr_multiply"
                                constantValues:values
                                         error:&error];
      if (!fn)
        throw std::runtime_error([[error localizedDescription] UTF8String]);
      pipeline = [m.device newComputePipelineStateWithFunction:fn error:&error];
      if (!pipeline)
        throw std::runtime_error([[error localizedDescription] UTF8String]);
      threads = std::min<size_t>(128, [pipeline maxTotalThreadsPerThreadgroup]);
      compile_ms = seconds(t) * 1000;
      t = Clock::now();
      auto put = [&](const void *data, size_t bytes) {
        id<MTLBuffer> b =
            [m.device newBufferWithLength:std::max<size_t>(bytes, 8)
                                  options:MTLResourceStorageModeShared];
        if (!b)
          throw std::runtime_error("Metal buffer allocation failed");
        if (data && bytes)
          std::memcpy([b contents], data, bytes);
        buffers.push_back(b);
      };
      put(p.row.data(), p.row.size() * 4);
      put(p.col.data(), p.col.size() * 4);
      put(p.a.data(), p.a.size() * 8);
      put(p.x.data(), p.x.size() * 8);
      put(nullptr, outputs * 8);
      put(p.perm.data(), p.perm.size() * 4);
      put(p.phase.data(), p.phase.size() * 8);
      upload_ms = seconds(t) * 1000;
    }
    Z *output() { return static_cast<Z *>([buffers[4] contents]); }
    double run(int repeats) {
      @autoreleasepool {
        id<MTLCommandBuffer> cmd = [m.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pipeline];
        for (size_t i = 0; i < buffers.size(); ++i)
          [enc setBuffer:buffers[i] offset:0 atIndex:i];
        [enc setBytes:&rows length:sizeof(rows) atIndex:7];
        for (int j = 0; j < repeats; ++j) {
          if (j)
            [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
          [enc dispatchThreads:MTLSizeMake(outputs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
        }
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
        if (cmd.status == MTLCommandBufferStatusError)
          throw std::runtime_error(
              [[cmd.error localizedDescription] UTF8String]);
        return ([cmd GPUEndTime] - [cmd GPUStartTime]) * 1000 / repeats;
      }
    }
  };

  struct Timing {
    double median_ms, min_ms, max_ms, device_ms = 0;
    int repeats;
  };
  Timing time_cpu(const std::function<void()> &f, int samples,
                  double target_ms) {
    f();
    f();
    auto t = Clock::now();
    f();
    double one = seconds(t) * 1000;
    int n = std::clamp(int(target_ms / std::max(one, 0.001)), 1, 128);
    std::vector<double> times;
    for (int s = 0; s < samples; ++s) {
      t = Clock::now();
      for (int i = 0; i < n; ++i)
        f();
      times.push_back(seconds(t) * 1000 / n);
    }
    std::sort(times.begin(), times.end());
    return {times[times.size() / 2], times.front(), times.back(), 0, n};
  }
  Timing time_gpu(GPUProblem &g, int samples, double target_ms, bool latency) {
    g.run(1);
    g.run(1);
    auto t = Clock::now();
    g.run(1);
    double one = seconds(t) * 1000;
    int n =
        latency ? 1 : std::clamp(int(target_ms / std::max(one, 0.001)), 1, 128);
    std::vector<double> times, device;
    for (int s = 0; s < samples; ++s) {
      t = Clock::now();
      double d = g.run(n);
      times.push_back(seconds(t) * 1000 / n);
      device.push_back(d);
    }
    std::sort(times.begin(), times.end());
    std::sort(device.begin(), device.end());
    return {times[times.size() / 2], times.front(), times.back(),
            device[device.size() / 2], n};
  }

  std::vector<uint32_t> parse_list(const std::string &s) {
    std::vector<uint32_t> out;
    std::stringstream in(s);
    std::string part;
    while (std::getline(in, part, ',')) {
      int n = std::stoi(part);
      if (n < 1 || n > 128)
        throw std::runtime_error("list values must be 1..128");
      out.push_back(uint32_t(n));
    }
    if (out.empty())
      throw std::runtime_error("empty list");
    return out;
  }

} // namespace MetalCloverEnv
