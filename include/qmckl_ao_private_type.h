#ifndef QMCKL_AO_HPT
#define QMCKL_AO_HPT

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

#include <omp.h>

#include "qmckl_blas_private_type.h"

/* Data structure                                                 :noexport: */

typedef struct qmckl_ao_basis_struct {
  int64_t shell_num;
  int64_t prim_num;
  int64_t ao_num;
  int64_t *restrict nucleus_index;
  int64_t *restrict nucleus_shell_num;
  int32_t *restrict shell_ang_mom;
  int64_t *restrict shell_prim_num;
  int64_t *restrict shell_prim_index;
  double *restrict shell_factor;
  double *restrict exponent;
  double *restrict coefficient;
  double *restrict prim_factor;
  double *restrict ao_factor;

  int64_t *restrict nucleus_prim_index;
  double *restrict coefficient_normalized;
  int32_t *restrict nucleus_max_ang_mom;
  double *restrict nucleus_range;
  double *restrict primitive_vgl;
  uint64_t primitive_vgl_date;
  double *restrict shell_vgl;
  uint64_t shell_vgl_date;
  double *restrict ao_vgl;
  uint64_t ao_vgl_date;
  double *restrict ao_value;
  uint64_t ao_value_date;

  int32_t uninitialized;
  bool provided;
  bool ao_cartesian;
  char type;
#ifdef HAVE_HPC
  /* HPC specific data structures */
  int32_t *restrict prim_num_per_nucleus;
  qmckl_tensor coef_per_nucleus;
  qmckl_matrix expo_per_nucleus;

#ifdef HAVE_DEVICE_POINTERS
  int32_t *restrict prim_num_per_nucleus_device;
#endif
#endif

#ifdef HAVE_DEVICE_POINTERS
  int64_t *nucleus_index_device;
  int64_t *nucleus_shell_num_device;
  int32_t *shell_ang_mom_device;
  int64_t *shell_prim_num_device;
  int64_t *shell_prim_index_device;
  double *shell_factor_device;
  double *exponent_device;
  double *coefficient_device;
  double *prim_factor_device;
  double *ao_factor_device;

  int64_t *nucleus_prim_index_device;
  double *coefficient_normalized_device;
  int32_t *nucleus_max_ang_mom_device;
  double *nucleus_range_device;
  double *primitive_vgl_device;
  double *shell_vgl_device;
  double *ao_vgl_device;
  double *ao_value_device;
#endif

} qmckl_ao_basis_struct;

#endif
