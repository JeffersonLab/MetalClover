#ifndef METALCLOVER_H
#define METALCLOVER_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
namespace MetalCloverEnv {
  struct Context;
  extern "C" {
#endif

#define METALCLOVER_ABI_VERSION 1
#if defined(__GNUC__)
#define METALCLOVER_API __attribute__((visibility("default")))
#else
#define METALCLOVER_API
#endif
/* All fields use interleaved FP64 {real,imag}. Lexicographic site order,
 * x fastest, then y,z,t; spinors: site,color,spin; links: site,mu,row,col.
 * Odd fields enumerate odd sites in this same order. One RHS per call.
 * Contexts own their gauge/hierarchy and are serial, not thread safe.
 */
#ifdef __cplusplus
  typedef MetalCloverEnv::Context metalclover_context;
#else
typedef struct metalclover_context metalclover_context;
#endif
  typedef enum {
    METALCLOVER_BICGSTAB = 0,
    METALCLOVER_MG = 1
  } metalclover_method;
  typedef struct {
    uint32_t abi_version, shape[4];
    double mass, xi0, nu, clover_r, clover_t;
    /* 1: un-BC'd stout links, apply temporal AP. 0: state links already BC'd.
     * The library never applies stout smearing to the supplied links. */
    int apply_temporal_antiperiodic, method, workers, restart;
    uint32_t blocks[4];
    int null_candidates, null_iterations, null_restart, null_warm, null_batch;
    double null_tolerance;
    int second_filter, bottom_iterations;
  } metalclover_options;
  typedef struct {
    int converged, outer_iterations, krylov_iterations;
    uint64_t fine_schur_applications, nominal_operator_flops;
    uint64_t fine_flops, coarse_flops, transfer_flops, local_vector_flops;
    double seconds, residual_seconds, relative_residual;
    /* Nominal dispatched operator work / complete solve time. Includes padded
     * coarse work; excludes Krylov vector ops and CPU FP64 work from numerator.
     * Denominator includes FP64 checks. This is not a hardware-counter rate. */
    double nominal_operator_gflops;
  } metalclover_solve_info;
  /* Optional independent CPU FP64 operator, PLUS sign, matching the solve
   * subset. Returns zero on success. Input/output each contain n_complex
   * interleaved complex numbers; access with memcpy or the caller's own type.
   */
  typedef int (*metalclover_apply64)(void *user, const void *input,
                                     void *output, size_t n_complex);
  METALCLOVER_API void
  metalclover_default_options(metalclover_options *options);
  METALCLOVER_API int metalclover_create(const metalclover_options *options,
                                         const void *stout_links64,
                                         metalclover_context **context);
  METALCLOVER_API void metalclover_destroy(metalclover_context *context);
  METALCLOVER_API const char *metalclover_last_error(void);
  METALCLOVER_API double
  metalclover_setup_seconds(const metalclover_context *context);
  /* Both solves use x as an initial guess and return info, without printing
   * solve results. max_iterations is a total Krylov budget over refinements.
   * A nonzero status means failure, including failure of the FP64 tolerance. */
  METALCLOVER_API int
  metalclover_solve_full(metalclover_context *context, const void *b64,
                         void *x64, double tolerance, int max_iterations,
                         metalclover_apply64 apply, void *user,
                         metalclover_solve_info *info);
  /* b64 is ALREADY reduced: S_oo x_o = b_o; no RHS reduction is repeated. */
  METALCLOVER_API int
  metalclover_solve_odd(metalclover_context *context, const void *b64,
                        void *x64, double tolerance, int max_iterations,
                        metalclover_apply64 apply, void *user,
                        metalclover_solve_info *info);
  /* Array calls reuse one context/hierarchy and one resident GPU RHS. Fields
   * are pointer arrays, not one contiguous multi-RHS allocation. Outputs must
   * be distinct and disjoint from all source fields. Each info corresponds to
   * the same array index. On failure, completed results remain available and
   * failed_rhs identifies the failed column; SIZE_MAX means all succeeded.
   * An empty array succeeds without touching fields or info. */
  METALCLOVER_API int metalclover_solve_full_array(
      metalclover_context *context, size_t rhs_count, const void *const *b64,
      void *const *x64, double tolerance, int max_iterations,
      metalclover_apply64 apply, void *user, metalclover_solve_info *info,
      size_t *failed_rhs);
  METALCLOVER_API int
  metalclover_solve_odd_array(metalclover_context *context, size_t rhs_count,
                              const void *const *b64, void *const *x64,
                              double tolerance, int max_iterations,
                              metalclover_apply64 apply, void *user,
                              metalclover_solve_info *info, size_t *failed_rhs);
  /* CPU FP64 helpers: b'_o=b_o-D_oe A_ee^-1 b_e;
   * x_e=A_ee^-1(b_e-D_eo x_o). They preserve full double precision. */
  METALCLOVER_API int metalclover_reduce_rhs(metalclover_context *context,
                                             const void *full_b64,
                                             void *odd_b64);
  METALCLOVER_API int metalclover_reconstruct(metalclover_context *context,
                                              const void *full_b64,
                                              const void *odd_x64,
                                              void *full_x64);
  METALCLOVER_API int metalclover_apply_full(metalclover_context *context,
                                             const void *x64, void *y64);
  METALCLOVER_API int metalclover_apply_odd(metalclover_context *context,
                                            const void *x64, void *y64);
#ifdef __cplusplus
  } // extern "C"
} // namespace MetalCloverEnv
#endif
#endif
