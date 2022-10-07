#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>


#include <qmckl.h>

#include "include/qmckl_context_private_type.h"
#include "include/qmckl_memory_private_type.h"
#include "include/qmckl_blas_private_type.h"
#include "include/qmckl_memory_private_func.h"
#include "include/qmckl_ao_private_type.h"
#include "include/qmckl_ao_private_func.h"



// qmckl_ao_openmp.c

qmckl_exit_code
qmckl_ao_polynomial_transp_vgl_hpc_omp_offload (const qmckl_context context,
                                    const double* restrict X,
                                    const double* restrict R,
                                    const int32_t lmax,
                                    int64_t* restrict n,
                                    const int64_t ldl,
                                    double* restrict const VGL,
                                    const int64_t ldv );

qmckl_exit_code
qmckl_compute_ao_vgl_gaussian_omp_offload (
                                           const qmckl_context context,
                                           const int64_t ao_num,
                                           const int64_t shell_num,
                                           const int32_t* restrict prim_num_per_nucleus,
                                           const int64_t point_num,
                                           const int64_t nucl_num,
                                           const double* restrict coord,
                                           const double* restrict nucl_coord,
                                           const int64_t* restrict nucleus_index,
                                           const int64_t* restrict nucleus_shell_num,
                                           const double* nucleus_range,
                                           const int32_t* restrict nucleus_max_ang_mom,
                                           const int32_t* restrict shell_ang_mom,
                                           const double* restrict ao_factor,
                                           const qmckl_matrix expo_per_nucleus,
                                           const qmckl_tensor coef_per_nucleus,
                                           double* restrict const ao_vgl );

qmckl_exit_code qmckl_provide_ao_vgl_omp_offload(qmckl_context context);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_omp_offload (qmckl_context context,
                                       double* const ao_vgl,
                                       const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace_omp_offload (qmckl_context context,
                                               double* const ao_vgl,
                                               const int64_t size_max);

// qmckl_ao_acc.c

qmckl_exit_code
qmckl_ao_polynomial_transp_vgl_hpc_acc_offload (const qmckl_context context,
                                    const double* restrict X,
                                    const double* restrict R,
                                    const int32_t lmax,
                                    int64_t* restrict n,
                                    const int64_t ldl,
                                    double* restrict const VGL,
                                    const int64_t ldv );

qmckl_exit_code
qmckl_compute_ao_vgl_gaussian_acc_offload (
                                           const qmckl_context context,
                                           const int64_t ao_num,
                                           const int64_t shell_num,
                                           const int32_t* restrict prim_num_per_nucleus,
                                           const int64_t point_num,
                                           const int64_t nucl_num,
                                           const double* restrict coord,
                                           const double* restrict nucl_coord,
                                           const int64_t* restrict nucleus_index,
                                           const int64_t* restrict nucleus_shell_num,
                                           const double* nucleus_range,
                                           const int32_t* restrict nucleus_max_ang_mom,
                                           const int32_t* restrict shell_ang_mom,
                                           const double* restrict ao_factor,
                                           const qmckl_matrix expo_per_nucleus,
                                           const qmckl_tensor coef_per_nucleus,
                                           double* restrict const ao_vgl );


qmckl_exit_code qmckl_provide_ao_vgl_acc_offload(qmckl_context context);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_acc_offload (qmckl_context context,
                                       double* const ao_vgl,
                                       const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace_acc_offload (qmckl_context context,
                                               double* const ao_vgl,
                                               const int64_t size_max);

// Device data structures

typedef struct qmckl_ao_basis_struct_device {
  int64_t            shell_num;
  int64_t            prim_num;
  int64_t            ao_num;
  int64_t * restrict nucleus_index;
  int64_t * restrict nucleus_shell_num;
  int32_t * restrict shell_ang_mom;
  int64_t * restrict shell_prim_num;
  int64_t * restrict shell_prim_index;
  double  * restrict shell_factor;
  double  * restrict exponent;
  double  * restrict coefficient;
  double  * restrict prim_factor;
  double  * restrict ao_factor;

  int64_t * restrict nucleus_prim_index;
  double  * restrict coefficient_normalized;
  int32_t * restrict nucleus_max_ang_mom;
  double  * restrict nucleus_range;
  double  * restrict primitive_vgl;
  uint64_t           primitive_vgl_date;
  double  * restrict shell_vgl;
  uint64_t           shell_vgl_date;
  double  * restrict ao_vgl;
  uint64_t           ao_vgl_date;
  double  * restrict ao_value;
  uint64_t           ao_value_date;

  int32_t            uninitialized;
  bool               provided;
  bool               ao_cartesian;
  char               type;

#ifdef HAVE_HPC
 int32_t*  restrict   prim_num_per_nucleus;
 qmckl_tensor         coef_per_nucleus;
 qmckl_matrix         expo_per_nucleus;

 int32_t*  restrict   prim_num_per_nucleus_device;
#endif

  int64_t * restrict nucleus_index_device;
  int64_t * restrict nucleus_shell_num_device;
  int32_t * restrict shell_ang_mom_device;
  int64_t * restrict shell_prim_num_device;
  int64_t * restrict shell_prim_index_device;
  double  * restrict shell_factor_device;
  double  * restrict exponent_device;
  double  * restrict coefficient_device;
  double  * restrict prim_factor_device;
  double  * restrict ao_factor_device;

  int64_t * restrict nucleus_prim_index_device;
  double  * restrict coefficient_normalized_device;
  int32_t * restrict nucleus_max_ang_mom_device;
  double  * restrict nucleus_range_device;
  double  * restrict primitive_vgl_device;
  double  * restrict shell_vgl_device;
  double  * restrict ao_vgl_device;
  double  * restrict ao_value_device;
} qmckl_ao_basis_struct_device;


// qmckl_ao_device.c

qmckl_exit_code
qmckl_ao_polynomial_transp_vgl_hpc_device (const qmckl_context context,
                                    const double* restrict X,
                                    const double* restrict R,
                                    const int32_t lmax,
                                    int64_t* restrict n,
                                    const int64_t ldl,
                                    double* restrict const VGL,
                                   const int64_t ldv );

qmckl_exit_code
qmckl_compute_ao_vgl_gaussian_device_pointers (
                                           const qmckl_context context,
                                           const int64_t ao_num,
                                           const int64_t shell_num,
                                           const int32_t* restrict prim_num_per_nucleus,
                                           const int64_t point_num,
                                           const int64_t nucl_num,
                                           const double* restrict coord,
                                           const double* restrict nucl_coord,
                                           const int64_t* restrict nucleus_index,
                                           const int64_t* restrict nucleus_shell_num,
                                           const double* nucleus_range,
                                           const int32_t* restrict nucleus_max_ang_mom,
                                           const int32_t* restrict shell_ang_mom,
                                           const double* restrict ao_factor,
                                           const qmckl_matrix expo_per_nucleus,
                                           const qmckl_tensor coef_per_nucleus,
                                           double* restrict const ao_vgl,

                                           int device_id );


qmckl_exit_code qmckl_provide_ao_basis_ao_vgl_device(qmckl_context context, int device_id);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_device (qmckl_context context,
                                  double* const ao_vgl,
                                  const int64_t size_max,
                                  int device_id);
