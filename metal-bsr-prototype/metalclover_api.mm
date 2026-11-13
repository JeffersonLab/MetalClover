#include "include/metalclover/metalclover.h"
#include "low_memory_clover.hpp"
#include "metalclover_shaders.hpp"
#include "mg_real_tune.hpp"
#include "near_null_svd.hpp"
#include <memory>
#include <mutex>

namespace MetalCloverEnv {
  namespace {
    thread_local std::string last_error;
    template <class F> int protect(F f) {
      @autoreleasepool {
        try {
          f();
          last_error.clear();
          return 0;
        } catch (const std::exception &e) {
          last_error = e.what();
          return 1;
        }
      }
    }
    template <class T> void release_storage(std::vector<T> &v) {
      std::vector<T>().swap(v);
    }
    std::vector<D> copy_field(const void *p, size_t n) {
      if (!p)
        throw std::runtime_error("null field pointer");
      std::vector<D> a(n);
      std::memcpy(a.data(), p, n * sizeof(D));
      return a;
    }
    using Coarse = GPUBlockCoarseSolve<BlockMatrix>;
    using CoarseFG = GPUGenericFGMRES<GPUBlockOperator, TunedGPUCycle>;
    using FineFG = GPUGenericFGMRES<GPUCompactSchur, TunedGPUCycle>;
  } // namespace

  struct Context {
    metalclover_options options;
    std::unique_ptr<MappedField> gauge;
    std::unique_ptr<LowMemoryClover> c;
    std::unique_ptr<MetalRun> metal;
    std::unique_ptr<GPUCompactSchur> op;
    std::unique_ptr<Transfer> first, second;
    std::unique_ptr<BlockMatrix> a1, a2;
    std::unique_ptr<TunedGPUCycle> pre, pre1;
    std::unique_ptr<Coarse> matrix1, bottom;
    std::unique_ptr<GPUBlockOperator> op1;
    std::unique_ptr<CoarseFG> coarse;
    std::unique_ptr<FineFG> mg;
    std::unique_ptr<GPUBiCGStab> bi;
    std::unique_ptr<CompactDeviceCorrection<FineFG>> mg_correction;
    std::unique_ptr<CompactDeviceCorrection<GPUBiCGStab>> bi_correction;

    struct BatchWorkspace {
      uint32_t rhs;
      GPUCompactSchur op;
      std::unique_ptr<TunedGPUCycle> pre, pre1;
      std::unique_ptr<Coarse> matrix1, bottom;
      std::unique_ptr<GPUBlockOperator> op1;
      std::unique_ptr<CoarseFG> coarse;
      std::unique_ptr<FineFG> mg;
      std::unique_ptr<GPUBiCGStab> bi;
      std::unique_ptr<CompactDeviceCorrection<FineFG>> mg_correction;
      std::unique_ptr<CompactDeviceCorrection<GPUBiCGStab>> bi_correction;
      size_t allocated_bytes = 0;
      BatchWorkspace(Context &context, uint32_t width);
    };
    std::unique_ptr<BatchWorkspace> batch;
    size_t batch_limit = 8;
    metalclover_batch_info batch_info{};
    double setup_seconds = 0;
    std::mutex mutex;
    uint32_t batch_width(size_t remaining) const;
    void solve_many(bool full, size_t count, const void *const *b,
                    void *const *x, double tolerance, int limit,
                    metalclover_apply64 apply, void *user,
                    metalclover_solve_info *info, size_t *failed);
    void solve_tile(bool full, size_t offset, uint32_t width,
                    const void *const *b, void *const *x, double tolerance,
                    int limit, metalclover_apply64 apply, void *user,
                    metalclover_solve_info *info, size_t *failed);

    explicit Context(const metalclover_options &o, const void *links)
        : options(o) {
      auto started = Clock::now();
      if (o.abi_version != METALCLOVER_ABI_VERSION || !links || o.workers < 1 ||
          o.workers > 64 ||
          (o.method != METALCLOVER_MG && o.method != METALCLOVER_BICGSTAB))
        throw std::runtime_error("invalid MetalClover options");
      production::Header h;
      std::copy_n(o.shape, 4, h.shape.begin());
      h.mass = o.mass;
      h.xi0 = o.xi0;
      h.nu = o.nu;
      h.clover_r = o.clover_r;
      h.clover_t = o.clover_t;
      h.time_count = 1;
      h.validate();
      for (auto extent : h.shape)
        if (extent % 2)
          throw std::runtime_error("even lattice extents required");
      // File backing lets the original FP64 gauge be evicted and reread safely.
      char path[] = "/private/tmp/metalclover-gauge-XXXXXX";
      int fd = mkstemp(path);
      if (fd < 0)
        throw std::runtime_error("create gauge backing file");
      size_t bytes = size_t(h.volume()) * 36 * sizeof(D);
      int err = ftruncate(fd, bytes);
      close(fd);
      if (err) {
        unlink(path);
        throw std::runtime_error("size gauge backing file");
      }
      try {
        gauge = std::make_unique<MappedField>(path, bytes);
      } catch (...) {
        unlink(path);
        throw;
      }
      unlink(path);
      std::memcpy(gauge->data, links, bytes);
      gauge->drop(0, 0, true);
      c = std::make_unique<LowMemoryClover>(*gauge, h, 0, o.workers,
                                            o.apply_temporal_antiperiodic != 0);
      metal = std::make_unique<MetalRun>(metalclover_shader);
      op = std::make_unique<GPUCompactSchur>(*metal, *c);
      release_storage(c->links32);
      release_storage(c->oddA32);
      release_storage(c->inverseEven32);
      c->drop_coefficients();
      if (o.method == METALCLOVER_BICGSTAB) {
        bi = std::make_unique<GPUBiCGStab>(*metal, c->parity.oe, 1, o.workers);
        bi->encode_override = [this](id<MTLComputeCommandEncoder> e,
                                     id<MTLBuffer> x,
                                     id<MTLBuffer> y) { op->encode(e, x, y); };
        bi_correction = std::make_unique<CompactDeviceCorrection<GPUBiCGStab>>(
            *c, *op, *bi, o.workers);
      } else {
        if (o.restart < 2 || o.restart > 12 || o.second_filter < 1 ||
            o.bottom_iterations < 1)
          throw std::runtime_error("invalid MG iteration options");
        LatticeShape blocks{}, coarse_shape{};
        for (int mu = 0; mu < 4; ++mu) {
          blocks[mu] = o.blocks[mu];
          if (!blocks[mu] || h.shape[mu] % blocks[mu] ||
              (h.shape[mu] / blocks[mu]) % 2)
            throw std::runtime_error("invalid MG aggregation geometry");
          coarse_shape[mu] = h.shape[mu] / blocks[mu];
        }
        first =
            std::make_unique<Transfer>(c->parity.oe.p.rows, 12, 8, h.shape,
                                       blocks, true, o.workers, &c->parity.odd);
        NearNullOptions null;
        null.candidates = o.null_candidates;
        null.iterations = o.null_iterations;
        null.restart = o.null_restart;
        null.warm = o.null_warm;
        null.batch = o.null_batch;
        null.tolerance = o.null_tolerance;
        std::ostringstream vector_log, spectrum;
        build_gpu_near_null_svd(*first, *op, null, vector_log, spectrum);
        c->drop_coefficients();
        a1 = std::make_unique<BlockMatrix>(16, first->aggregates, o.workers);
        assemble_schur(*a1, *c, *first);
        c->drop_coefficients();
        second = std::make_unique<Transfer>(a1->rows, 16, 8, coarse_shape,
                                            hypercube(2), false, o.workers);
        second->build([&](const D *x, D *y) { a1->apply(x, y, 8); },
                      [&](const D *x, D *y) { a1->apply(x, y, 8, true); },
                      o.second_filter);
        a2 = std::make_unique<BlockMatrix>(16, second->aggregates, o.workers);
        project_matrix(*a2, *a1, *second);
        first->release64();
        second->release64();
        a1->release64();
        a2->release64();
        pre = std::make_unique<TunedGPUCycle>(*first, *metal, 1,
                                              std::vector<Cx<float>>{}, 2);
        pre->inv = op->Aodd;
        release_storage(first->basis32);
        pre->product = [this](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                              id<MTLBuffer> y) { op->encode(e, x, y); };
        matrix1 = std::make_unique<Coarse>(*a1, *metal, 1, 1, .1);
        bottom =
            std::make_unique<Coarse>(*a2, *metal, 1, o.bottom_iterations, .1);
        op1 = std::make_unique<GPUBlockOperator>(*matrix1);
        pre1 = std::make_unique<TunedGPUCycle>(*second, *metal, 1, a1->di32);
        release_storage(second->basis32);
        pre1->product = [this](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                               id<MTLBuffer> y) { op1->encode(e, x, y); };
        pre1->coarse = [this](id<MTLComputeCommandEncoder> e, id<MTLBuffer> x,
                              id<MTLBuffer> y) { bottom->encode(e, x, y); };
        pre1->check = [this]() { bottom->check(); };
        coarse = std::make_unique<CoarseFG>(a1->n, 1, *op1, *pre1, 8);
        coarse->allow_inexact = true;
        pre->blocking_coarse = [this](id<MTLBuffer> x, id<MTLBuffer> y) {
          coarse->run_device(x, .1, 32);
          pre->copy(coarse->delta, y, a1->n);
        };
        release_storage(a1->a32);
        release_storage(a1->di32);
        release_storage(a2->a32);
        release_storage(a2->di32);
        mg =
            std::make_unique<FineFG>(first->rows * 12, 1, *op, *pre, o.restart);
        mg->blocking_pre = [this](id<MTLBuffer> x, id<MTLBuffer> y) {
          pre->apply_device(x, y);
        };
        mg_correction = std::make_unique<CompactDeviceCorrection<FineFG>>(
            *c, *op, *mg, o.workers);
      }
      setup_seconds = seconds(started);
    }

    void local_even(const D *x, D *y) {
      chunks(c->parity.even.size(), options.workers,
             [&](int, uint32_t begin, uint32_t end) {
               for (uint32_t site = begin; site < end; ++site) {
                 D inv[72];
                 c->even_inverse_blocks(site, inv);
                 for (uint32_t ch = 0; ch < 2; ++ch)
                   for (uint32_t i = 0; i < 6; ++i) {
                     CD sum = 0;
                     for (uint32_t j = 0; j < 6; ++j) {
                       auto a = inv[ch * 36 + i * 6 + j],
                            v = x[size_t(site) * 12 + (j % 3) * 4 + ch * 2 +
                                  j / 3];
                       sum += CD(a.re, a.im) * CD(v.re, v.im);
                     }
                     y[size_t(site) * 12 + (i % 3) * 4 + ch * 2 + i / 3] = {
                         sum.real(), sum.imag()};
                   }
               }
             });
    }
    void apply_odd(const D *x, D *y) {
      size_t n = size_t(c->parity.odd.size()) * 12;
      std::vector<D> temp(n), even(n);
      c->apply(x, temp.data(), options.workers, 0, false);
      local_even(temp.data(), even.data());
      c->apply(even.data(), y, options.workers, 1, false);
      chunks(c->parity.odd.size(), options.workers,
             [&](int, uint32_t begin, uint32_t end) {
               for (uint32_t site = begin; site < end; ++site)
                 for (uint32_t ch = 0; ch < 2; ++ch)
                   for (uint32_t i = 0; i < 6; ++i) {
                     CD sum = 0;
                     for (uint32_t j = 0; j < 6; ++j) {
                       auto a = c->onsite_value(true, site, ch, i, j);
                       auto v =
                           x[size_t(site) * 12 + (j % 3) * 4 + ch * 2 + j / 3];
                       sum += CD(a.re, a.im) * CD(v.re, v.im);
                     }
                     size_t q =
                         size_t(site) * 12 + (i % 3) * 4 + ch * 2 + i / 3;
                     y[q] = {sum.real() - y[q].re, sum.imag() - y[q].im};
                   }
             });
    }
    void reduce(const D *b, D *out) {
      size_t n = size_t(c->parity.odd.size()) * 12;
      std::vector<D> even(n), inverse(n);
      for (size_t i = 0; i < c->parity.even.size(); ++i)
        std::copy_n(b + size_t(c->parity.even[i]) * 12, 12,
                    even.data() + i * 12);
      local_even(even.data(), inverse.data());
      c->apply(inverse.data(), out, options.workers, 1, false);
      for (size_t i = 0; i < c->parity.odd.size(); ++i)
        for (size_t q = 0; q < 12; ++q)
          out[i * 12 + q] =
              minus_z(b[size_t(c->parity.odd[i]) * 12 + q], out[i * 12 + q]);
    }
    void reconstruct(const D *b, const D *odd, D *out) {
      size_t n = size_t(c->parity.odd.size()) * 12;
      std::vector<D> even(n), inverse(n);
      c->apply(odd, even.data(), options.workers, 0, false);
      for (size_t i = 0; i < c->parity.even.size(); ++i)
        for (size_t q = 0; q < 12; ++q)
          even[i * 12 + q] =
              minus_z(b[size_t(c->parity.even[i]) * 12 + q], even[i * 12 + q]);
      local_even(even.data(), inverse.data());
      for (size_t i = 0; i < c->parity.even.size(); ++i) {
        std::copy_n(inverse.data() + i * 12, 12,
                    out + size_t(c->parity.even[i]) * 12);
        std::copy_n(odd + i * 12, 12, out + size_t(c->parity.odd[i]) * 12);
      }
    }
    void solve(bool full, const void *bptr, void *xptr, double tolerance,
               int limit, metalclover_apply64 apply, void *user,
               metalclover_solve_info &info, bool acquire_lock = true) {
      std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
      if (acquire_lock && !lock.try_lock())
        throw std::runtime_error("concurrent use of solver context");
      if (!(tolerance > 0 && tolerance < 1) || limit < 1)
        throw std::runtime_error("invalid solve tolerance or iteration limit");
      auto start = Clock::now();
      size_t n = size_t(c->w.a.p.rows) * 12 / (full ? 1 : 2);
      auto b = copy_field(bptr, n), x = copy_field(xptr, n);
      std::vector<D> r(n);
      auto before = metal->work;
      auto applications_before = op->applications;
      auto norm = [](const std::vector<D> &v) {
        double sum = 0;
        for (auto z : v)
          sum += z.re * z.re + z.im * z.im;
        return std::sqrt(sum);
      };
      double bnorm = norm(b);
      if (bnorm == 0)
        bnorm = 1;
      auto check = [&]() {
        auto t = Clock::now();
        if (apply) {
          if (apply(user, x.data(), r.data(), n))
            throw std::runtime_error("external FP64 operator callback failed");
        } else if (full)
          c->apply(x.data(), r.data(), options.workers, -1, true);
        else
          apply_odd(x.data(), r.data());
        for (size_t j = 0; j < n; ++j)
          r[j] = minus_z(b[j], r[j]);
        info.relative_residual = norm(r) / bnorm;
        info.residual_seconds += seconds(t);
        if (!std::isfinite(info.relative_residual))
          throw std::runtime_error("nonfinite FP64 residual");
      };
      check();
      while (info.relative_residual > tolerance && info.outer_iterations < 16 &&
             info.krylov_iterations < limit) {
        int budget = limit - info.krylov_iterations;
        auto correction = [&](auto &corr) {
          InnerInfo inner;
          const Cx<float> *dx;
          // Normalize before FP32 conversion; otherwise independent RHS scales
          // can underflow FP32 norm reductions near a 1e-12 relative target.
          double scale = norm(r);
          auto *p = static_cast<Cx<float> *>(
              [(full ? corr.input : corr.reduced) contents]);
          for (size_t j = 0; j < n; ++j)
            p[j] = {float(r[j].re / scale), float(r[j].im / scale)};
          if (full) {
            inner = corr.run_prepared(.002, budget);
            dx = corr.correction();
          } else {
            inner = corr.solve_reduced(.002, budget);
            dx = corr.odd_correction();
          }
          for (size_t j = 0; j < n; ++j) {
            x[j].re += scale * dx[j].re;
            x[j].im += scale * dx[j].im;
          }
          info.krylov_iterations += inner.max_iterations;
        };
        if (mg_correction)
          correction(*mg_correction);
        else
          correction(*bi_correction);
        ++info.outer_iterations;
        check();
      }
      std::memcpy(xptr, x.data(), n * sizeof(D));
      info.seconds = seconds(start);
      auto work = metal->work - before;
      info.fine_schur_applications = op->applications - applications_before;
      info.nominal_operator_flops = work.total();
      info.fine_flops = work.fine;
      info.coarse_flops = work.coarse;
      info.transfer_flops = work.transfer;
      info.local_vector_flops = work.local_vector;
      info.nominal_operator_gflops = work.gflops(info.seconds);
      info.converged = info.relative_residual <= tolerance;
      if (!info.converged)
        throw std::runtime_error("FP64 residual tolerance not attained");
    }
  };

  Context::BatchWorkspace::BatchWorkspace(Context &ctx, uint32_t width)
      : rhs(width), op(*ctx.op, width) {
    // Coefficients and transfer weights are shared with the scalar graph.
    // Only RHS-dependent vectors and function specializations are new.
    if (ctx.options.method == METALCLOVER_BICGSTAB) {
      System system;
      system.p.rows = op.eo.rows;
      system.p.rhs = rhs;
      system.p.colors = 3;
      system.p.spins = 4;
      bi = std::make_unique<GPUBiCGStab>(*ctx.metal, system, 1,
                                         ctx.options.workers);
      bi->allow_inexact = true;
      bi->encode_override = [this](id<MTLComputeCommandEncoder> enc,
                                   id<MTLBuffer> x,
                                   id<MTLBuffer> y) { op.encode(enc, x, y); };
      bi_correction = std::make_unique<CompactDeviceCorrection<GPUBiCGStab>>(
          *ctx.c, op, *bi, ctx.options.workers);
    } else {
      std::vector<Cx<float>> empty;
      pre = std::make_unique<TunedGPUCycle>(*ctx.first, *ctx.metal, rhs, empty,
                                            2, ctx.pre.get());
      pre->product = [this](id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
                            id<MTLBuffer> y) { op.encode(enc, x, y); };
      matrix1 = std::make_unique<Coarse>(*ctx.a1, *ctx.metal, rhs, 1, .1,
                                         ctx.matrix1.get());
      bottom = std::make_unique<Coarse>(*ctx.a2, *ctx.metal, rhs,
                                        ctx.options.bottom_iterations, .1,
                                        ctx.bottom.get());
      op1 = std::make_unique<GPUBlockOperator>(*matrix1);
      pre1 = std::make_unique<TunedGPUCycle>(*ctx.second, *ctx.metal, rhs,
                                             empty, 0, ctx.pre1.get());
      pre1->product = [this](id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
                             id<MTLBuffer> y) { op1->encode(enc, x, y); };
      pre1->coarse = [this](id<MTLComputeCommandEncoder> enc, id<MTLBuffer> x,
                            id<MTLBuffer> y) { bottom->encode(enc, x, y); };
      pre1->check = [this]() { bottom->check(); };
      coarse = std::make_unique<CoarseFG>(ctx.a1->n, rhs, *op1, *pre1, 8);
      coarse->allow_inexact = true;
      pre->blocking_coarse = [this, &ctx](id<MTLBuffer> x, id<MTLBuffer> y) {
        coarse->run_device(x, .1, 32);
        pre->copy(coarse->delta, y, size_t(ctx.a1->n) * rhs);
      };
      mg = std::make_unique<FineFG>(ctx.first->rows * 12, rhs, op, *pre,
                                    ctx.options.restart);
      mg->allow_inexact = true;
      mg->blocking_pre = [this](id<MTLBuffer> x, id<MTLBuffer> y) {
        pre->apply_device(x, y);
      };
      mg_correction = std::make_unique<CompactDeviceCorrection<FineFG>>(
          *ctx.c, op, *mg, ctx.options.workers);
    }
  }

  uint32_t Context::batch_width(size_t remaining) const {
    size_t width = std::min(remaining, batch_limit);
    size_t components = size_t(c->w.a.p.rows) * 12;
    width = std::min(width, size_t(UINT32_MAX) / components);
    width = std::max(size_t(1), width);
    if (batch && batch->rhs == width)
      return uint32_t(width);
    // Keep the old large-volume scalar path viable: batching must fit in the
    // device's recommended budget in addition to the retained scalar graph.
    double per_rhs = double(components / 2) * 8;
    per_rhs *= options.method == METALCLOVER_MG ? 32 : 12;
    if (options.method == METALCLOVER_MG)
      per_rhs += (36. * a1->n + 10. * a2->n) * 8;
    per_rhs *= 1.08; // norm/Arnoldi partials and small allocation padding
    double current = double(metal->device.currentAllocatedSize);
    if (batch)
      current -= batch->allocated_bytes;
    double available = double(metal->device.recommendedMaxWorkingSetSize) -
                       current - 256. * 1024 * 1024;
    size_t fitting = available > 0 ? size_t(available / per_rhs) : 0;
    return uint32_t(std::max(size_t(1), std::min(width, fitting)));
  }

  void Context::solve_tile(bool full, size_t offset, uint32_t width,
                           const void *const *bptr, void *const *xptr,
                           double tolerance, int limit,
                           metalclover_apply64 apply, void *user,
                           metalclover_solve_info *info, size_t *failed) {
    if (width == 1) {
      if (failed)
        *failed = offset;
      auto before = op->applications;
      solve(full, bptr[offset], xptr[offset], tolerance, limit, apply, user,
            info[offset], false);
      batch_info.fine_schur_dispatches += op->applications - before;
      return;
    }
    if (!batch || batch->rhs != width) {
      auto started = Clock::now();
      batch.reset();
      size_t before = metal->device.currentAllocatedSize;
      batch = std::make_unique<BatchWorkspace>(*this, width);
      batch->allocated_bytes = metal->device.currentAllocatedSize - before;
      batch_info.workspace_seconds += seconds(started);
    }
    struct Field {
      std::vector<D> b, x, residual;
      double bnorm = 1, scale = 1;
    };
    size_t n = size_t(c->w.a.p.rows) * 12 / (full ? 1 : 2);
    std::vector<Field> fields(width);
    auto norm = [](const std::vector<D> &v) {
      double sum = 0;
      for (auto z : v)
        sum += z.re * z.re + z.im * z.im;
      return std::sqrt(sum);
    };
    auto publish = [&](uint32_t r) {
      std::memcpy(xptr[offset + r], fields[r].x.data(), n * sizeof(D));
    };
    auto check_field = [&](uint32_t r) {
      if (failed)
        *failed = offset + r;
      auto started = Clock::now();
      auto &f = fields[r];
      auto &out = info[offset + r];
      if (apply) {
        if (apply(user, f.x.data(), f.residual.data(), n))
          throw std::runtime_error(
              "external FP64 batch operator callback failed");
      } else if (full)
        c->apply(f.x.data(), f.residual.data(), options.workers, -1, true);
      else
        apply_odd(f.x.data(), f.residual.data());
      for (size_t j = 0; j < n; ++j)
        f.residual[j] = minus_z(f.b[j], f.residual[j]);
      f.scale = norm(f.residual);
      out.relative_residual = f.scale / f.bnorm;
      out.residual_seconds += seconds(started);
      if (!std::isfinite(out.relative_residual))
        throw std::runtime_error("nonfinite full FP64 batch residual");
      out.converged = out.relative_residual <= tolerance;
      if (out.converged)
        publish(r);
    };
    for (uint32_t r = 0; r < width; ++r) {
      auto &f = fields[r];
      f.b = copy_field(bptr[offset + r], n);
      f.x = copy_field(xptr[offset + r], n);
      f.residual.resize(n);
      f.bnorm = norm(f.b);
      if (f.bnorm == 0)
        f.bnorm = 1;
      check_field(r);
    }
    for (;;) {
      std::vector<uint32_t> active;
      int budget = limit;
      for (uint32_t r = 0; r < width; ++r) {
        auto &out = info[offset + r];
        if (out.converged)
          continue;
        if (failed)
          *failed = offset + r;
        if (out.outer_iterations >= 16 || out.krylov_iterations >= limit) {
          publish(r);
          throw std::runtime_error(
              "batch FP64 residual tolerance not attained within budget");
        }
        active.push_back(r);
        budget = std::min(budget, limit - out.krylov_iterations);
      }
      if (active.empty())
        break;
      if (failed)
        *failed = offset + active.front();
      auto before = metal->work;
      uint64_t applications_before = batch->op.applications;
      InnerInfo inner;
      const Cx<float> *dx = nullptr;
      auto correct = [&](auto &corr) {
        auto *input = static_cast<Cx<float> *>(
            [(full ? corr.input : corr.reduced) contents]);
        chunks(uint32_t(n), options.workers,
               [&](int, uint32_t begin, uint32_t end) {
                 for (size_t q = begin; q < end; ++q)
                   for (uint32_t r = 0; r < width; ++r) {
                     auto z = info[offset + r].converged
                                  ? D{}
                                  : fields[r].residual[q];
                     double scale =
                         info[offset + r].converged ? 1 : fields[r].scale;
                     input[q * width + r] = {float(z.re / scale),
                                             float(z.im / scale)};
                   }
               });
        inner = full ? corr.run_prepared(.002, budget)
                     : corr.solve_reduced(.002, budget);
        dx = full ? corr.correction() : corr.odd_correction();
      };
      if (batch->mg_correction)
        correct(*batch->mg_correction);
      else
        correct(*batch->bi_correction);
      if (inner.rhs_iterations.size() != width)
        throw std::runtime_error("missing per-RHS batch iteration counts");
      auto work = metal->work - before;
      uint64_t dispatches = batch->op.applications - applications_before;
      batch_info.fine_schur_dispatches += dispatches;
      // Attribute all nominal dispatched work, including inactive lanes, over
      // the columns corrected in this outer step. Preserve exact-guess zeros.
      auto share = [&](uint64_t total, size_t index) {
        return total / active.size() + uint64_t(index < total % active.size());
      };
      for (size_t k = 0; k < active.size(); ++k) {
        uint32_t r = active[k];
        auto &out = info[offset + r];
        out.fine_flops += share(work.fine, k);
        out.coarse_flops += share(work.coarse, k);
        out.transfer_flops += share(work.transfer, k);
        out.local_vector_flops += share(work.local_vector, k);
        out.nominal_operator_flops = out.fine_flops + out.coarse_flops +
                                     out.transfer_flops +
                                     out.local_vector_flops;
        out.fine_schur_applications += dispatches;
        out.krylov_iterations += inner.rhs_iterations[r];
        ++out.outer_iterations;
        auto &f = fields[r];
        for (size_t q = 0; q < n; ++q) {
          f.x[q].re += f.scale * dx[q * width + r].re;
          f.x[q].im += f.scale * dx[q * width + r].im;
        }
        check_field(r);
      }
    }
  }

  void Context::solve_many(bool full, size_t count, const void *const *b,
                           void *const *x, double tolerance, int limit,
                           metalclover_apply64 apply, void *user,
                           metalclover_solve_info *info, size_t *failed) {
    std::unique_lock<std::mutex> lock(mutex, std::try_to_lock);
    if (!lock.owns_lock())
      throw std::runtime_error("concurrent use of solver context");
    batch_info = {};
    batch_info.rhs_count = count;
    batch_info.configured_batch_rhs = batch_limit;
    auto started = Clock::now();
    auto before = metal->work;
    auto finish = [&]() {
      batch_info.seconds = seconds(started);
      batch_info.nominal_operator_flops = (metal->work - before).total();
      batch_info.nominal_operator_gflops =
          (metal->work - before).gflops(batch_info.seconds);
      for (size_t i = 0; i < count; ++i) {
        info[i].seconds = batch_info.seconds / count;
        info[i].nominal_operator_gflops =
            info[i].seconds > 0 ? double(info[i].nominal_operator_flops) /
                                      (info[i].seconds * 1e9)
                                : 0;
        batch_info.residual_seconds += info[i].residual_seconds;
      }
    };
    try {
      if (count && (!(tolerance > 0 && tolerance < 1) || limit < 1))
        throw std::runtime_error("invalid batch solve stopping parameters");
      for (size_t i = 0; i < count;) {
        uint32_t width = batch_width(count - i);
        ++batch_info.tiles;
        batch_info.maximum_resident_rhs =
            std::max(batch_info.maximum_resident_rhs, size_t(width));
        solve_tile(full, i, width, b, x, tolerance, limit, apply, user, info,
                   failed);
        i += width;
      }
      if (failed)
        *failed = SIZE_MAX;
    } catch (...) {
      finish();
      throw;
    }
    finish();
  }

  int solve_array(Context *context, bool full, size_t count,
                  const void *const *b, void *const *x, double tolerance,
                  int limit, metalclover_apply64 apply, void *user,
                  metalclover_solve_info *info, size_t *failed) {
    if (failed)
      *failed = 0;
    return protect([&] {
      if (!context || (count && (!b || !x || !info)))
        throw std::runtime_error("null array solve argument");
      size_t bytes =
          size_t(context->c->w.a.p.rows) * 12 * sizeof(D) / (full ? 1 : 2);
      auto overlap = [&](const void *a, const void *b) {
        auto left = reinterpret_cast<uintptr_t>(a),
             right = reinterpret_cast<uintptr_t>(b);
        if (left > UINTPTR_MAX - bytes || right > UINTPTR_MAX - bytes)
          throw std::runtime_error("invalid field address range");
        return left < right + bytes && right < left + bytes;
      };
      for (size_t i = 0; i < count; ++i)
        info[i] = {};
      for (size_t i = 0; i < count; ++i) {
        if (failed)
          *failed = i;
        if (!b[i] || !x[i])
          throw std::runtime_error("null array field");
        for (size_t j = 0; j < count; ++j)
          if (overlap(x[i], b[j]) || (i != j && overlap(x[i], x[j])))
            throw std::runtime_error("array output fields must be distinct and "
                                     "disjoint from sources");
      }
      context->solve_many(full, count, b, x, tolerance, limit, apply, user,
                          info, failed);
    });
  }

  extern "C" {
  int metalclover_set_rhs_batch(metalclover_context *ctx, size_t cap) {
    return protect([&] {
      if (!ctx || cap < 1 || cap > 64)
        throw std::runtime_error("RHS batch cap must be between 1 and 64");
      std::unique_lock<std::mutex> lock(ctx->mutex, std::try_to_lock);
      if (!lock.owns_lock())
        throw std::runtime_error("concurrent use of solver context");
      ctx->batch_limit = cap;
    });
  }
  int metalclover_get_batch_info(const metalclover_context *ctx,
                                 metalclover_batch_info *info) {
    return protect([&] {
      if (!ctx || !info)
        throw std::runtime_error("null batch info argument");
      *info = ctx->batch_info;
    });
  }
  void metalclover_default_options(metalclover_options *o) {
    if (!o)
      return;
    *o = {};
    o->abi_version = METALCLOVER_ABI_VERSION;
    o->xi0 = o->nu = o->clover_r = o->clover_t = 1;
    o->apply_temporal_antiperiodic = 1;
    o->method = METALCLOVER_MG;
    o->workers = 8;
    o->restart = 8;
    for (auto &b : o->blocks)
      b = 4;
    o->null_candidates = 24;
    o->null_iterations = 128;
    o->null_restart = 8;
    o->null_warm = 32;
    o->null_batch = 4;
    o->null_tolerance = .001;
    o->second_filter = 64;
    o->bottom_iterations = 8;
  }
  int metalclover_create(const metalclover_options *o, const void *g,
                         metalclover_context **ctx) {
    if (ctx)
      *ctx = nullptr;
    return protect([&] {
      if (!o || !ctx)
        throw std::runtime_error("null creation argument");
      auto p = std::make_unique<metalclover_context>(*o, g);
      *ctx = p.release();
    });
  }
  void metalclover_destroy(metalclover_context *ctx) { delete ctx; }
  const char *metalclover_last_error(void) { return last_error.c_str(); }
  double metalclover_setup_seconds(const metalclover_context *ctx) {
    return ctx ? ctx->setup_seconds : 0;
  }
  int metalclover_solve_full(metalclover_context *ctx, const void *b, void *x,
                             double tol, int limit, metalclover_apply64 apply,
                             void *user, metalclover_solve_info *info) {
    if (info)
      *info = {};
    return protect([&] {
      if (!ctx || !info)
        throw std::runtime_error("null solve argument");
      ctx->solve(true, b, x, tol, limit, apply, user, *info);
    });
  }
  int metalclover_solve_odd(metalclover_context *ctx, const void *b, void *x,
                            double tol, int limit, metalclover_apply64 apply,
                            void *user, metalclover_solve_info *info) {
    if (info)
      *info = {};
    return protect([&] {
      if (!ctx || !info)
        throw std::runtime_error("null solve argument");
      ctx->solve(false, b, x, tol, limit, apply, user, *info);
    });
  }
  int metalclover_solve_full_array(metalclover_context *context, size_t count,
                                   const void *const *b, void *const *x,
                                   double tolerance, int limit,
                                   metalclover_apply64 apply, void *user,
                                   metalclover_solve_info *info,
                                   size_t *failed) {
    return solve_array(context, true, count, b, x, tolerance, limit, apply,
                       user, info, failed);
  }
  int metalclover_solve_odd_array(metalclover_context *context, size_t count,
                                  const void *const *b, void *const *x,
                                  double tolerance, int limit,
                                  metalclover_apply64 apply, void *user,
                                  metalclover_solve_info *info,
                                  size_t *failed) {
    return solve_array(context, false, count, b, x, tolerance, limit, apply,
                       user, info, failed);
  }
  int metalclover_reduce_rhs(metalclover_context *ctx, const void *b,
                             void *out) {
    return protect([&] {
      if (!ctx || !out)
        throw std::runtime_error("null reduction argument");
      auto v = copy_field(b, size_t(ctx->c->w.a.p.rows) * 12);
      std::vector<D> r(v.size() / 2);
      ctx->reduce(v.data(), r.data());
      std::memcpy(out, r.data(), r.size() * sizeof(D));
    });
  }
  int metalclover_reconstruct(metalclover_context *ctx, const void *b,
                              const void *odd, void *out) {
    return protect([&] {
      if (!ctx || !out)
        throw std::runtime_error("null reconstruction argument");
      auto v = copy_field(b, size_t(ctx->c->w.a.p.rows) * 12),
           o = copy_field(odd, v.size() / 2);
      std::vector<D> x(v.size());
      ctx->reconstruct(v.data(), o.data(), x.data());
      std::memcpy(out, x.data(), x.size() * sizeof(D));
    });
  }
  int metalclover_apply_full(metalclover_context *ctx, const void *x,
                             void *out) {
    return protect([&] {
      if (!ctx || !out)
        throw std::runtime_error("null apply argument");
      auto v = copy_field(x, size_t(ctx->c->w.a.p.rows) * 12);
      std::vector<D> y(v.size());
      ctx->c->apply(v.data(), y.data(), ctx->options.workers, -1, true);
      std::memcpy(out, y.data(), y.size() * sizeof(D));
    });
  }
  int metalclover_apply_odd(metalclover_context *ctx, const void *x,
                            void *out) {
    return protect([&] {
      if (!ctx || !out)
        throw std::runtime_error("null apply argument");
      auto v = copy_field(x, size_t(ctx->c->w.a.p.rows) * 6);
      std::vector<D> y(v.size());
      ctx->apply_odd(v.data(), y.data());
      std::memcpy(out, y.data(), y.size() * sizeof(D));
    });
  }
  }

} // namespace MetalCloverEnv
