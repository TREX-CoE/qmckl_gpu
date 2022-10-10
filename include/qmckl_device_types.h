// This file is where custom GPU/device datatypes are defined

#pragma once

#include <stdint.h>
#include <qmckl.h>

#include "qmckl_numprec_private_type.h"


//**********
// BLAS types (matrices, ...)
//**********

// Vector
typedef struct qmckl_vector_device {
  double* restrict data;
  double* restrict data_device;

  int64_t size;
} qmckl_vector_device;

// Matrix
typedef struct qmckl_matrix_device {
  double* restrict data;
  double* restrict data_device;

  int64_t size[2];
} qmckl_matrix_device;

// Tensor
#define QMCKL_TENSOR_ORDER_MAX_DEVICE 16
typedef struct qmckl_tensor_device {
  double* restrict data;
  int64_t order;
  int64_t size[QMCKL_TENSOR_ORDER_MAX_DEVICE];

  double* restrict data_device;
} qmckl_tensor_device;


//**********
// Context subtypes
//**********

// Point
typedef struct qmckl_point_struct_device {
  int64_t             num;
  uint64_t            date;
  qmckl_matrix_device coord;
} qmckl_point_struct_device;

// Nucleus
typedef struct qmckl_nucleus_struct_device {
  int64_t      num;
  int64_t      repulsion_date;
  int64_t      nn_distance_date;
  int64_t      nn_distance_rescaled_date;
  int64_t      coord_date;
  qmckl_vector_device charge;
  qmckl_matrix_device coord;
  qmckl_matrix_device nn_distance;
  qmckl_matrix_device nn_distance_rescaled;
  double       repulsion;
  double       rescale_factor_kappa;
  int32_t      uninitialized;
  bool         provided;
} qmckl_nucleus_struct_device;

// Electron
typedef struct qmckl_walker_struct_device {
  int64_t num;
  qmckl_point_struct_device point;
} qmckl_walker_device;

typedef struct qmckl_electron_struct_device {
  int64_t        num;
  int64_t        up_num;
  int64_t        down_num;
  qmckl_walker_device   walker;
  qmckl_walker_device   walker_old;
  double         rescale_factor_kappa_ee;
  double         rescale_factor_kappa_en;
  uint64_t        ee_distance_date;
  uint64_t        en_distance_date;
  uint64_t        ee_potential_date;
  uint64_t        en_potential_date;
  uint64_t        ee_distance_rescaled_date;
  uint64_t        ee_distance_rescaled_deriv_e_date;
  uint64_t        en_distance_rescaled_date;
  uint64_t        en_distance_rescaled_deriv_e_date;
  double*        ee_distance;
  double*        en_distance;
  double*        ee_potential;
  double*        en_potential;
  double*        ee_distance_rescaled;
  double*        ee_distance_rescaled_deriv_e;
  double*        en_distance_rescaled;
  double*        en_distance_rescaled_deriv_e;
  int32_t        uninitialized;
  bool           provided;

  double*        ee_distance_device;
  double*        en_distance_device;
  double*        ee_potential_device;
  double*        en_potential_device;
  double*        ee_distance_rescaled_device;
  double*        ee_distance_rescaled_deriv_e_device;
  double*        en_distance_rescaled_device;
  double*        en_distance_rescaled_deriv_e_device;
} qmckl_electron_struct_device;

// AO
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

  int32_t*  restrict   prim_num_per_nucleus;
  int32_t*  restrict   prim_num_per_nucleus_device;

  qmckl_tensor_device  coef_per_nucleus;
  qmckl_matrix_device  expo_per_nucleus;

  int64_t * nucleus_index_device;
  int64_t * nucleus_shell_num_device;
  int32_t * shell_ang_mom_device;
  int64_t * shell_prim_num_device;
  int64_t * shell_prim_index_device;
  double  * shell_factor_device;
  double  * exponent_device;
  double  * coefficient_device;
  double  * prim_factor_device;
  double  * ao_factor_device;

  int64_t * nucleus_prim_index_device;
  double  * coefficient_normalized_device;
  int32_t * nucleus_max_ang_mom_device;
  double  * nucleus_range_device;
  double  * primitive_vgl_device;
  double  * shell_vgl_device;
  double  * ao_vgl_device;
  double  * ao_value_device;

} qmckl_ao_basis_struct_device;


//**********
// Main context type
//**********

typedef struct qmckl_context_struct_device {
  /* -- State of the library -- */

  /* Validity checking */
  uint64_t                 tag;

  /* Numerical precision */
  qmckl_numprec_struct     numprec;

  /* Thread lock */
  int                      lock_count;
  pthread_mutex_t          mutex;

  /* Error handling */
  qmckl_error_struct       error;

  /* Memory allocation */
  qmckl_memory_struct      memory;

  /* Current date */
  uint64_t                 date;

  /* Points */
  qmckl_point_struct_device         point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct_device       nucleus;
  qmckl_electron_struct_device      electron;
  qmckl_ao_basis_struct_device      ao_basis;
  qmckl_mo_basis_struct             mo_basis; // TODO
  qmckl_jastrow_struct              jastrow; // TODO
  qmckl_determinant_struct          det;
  qmckl_local_energy_struct         local_energy;

} qmckl_context_struct_device;


typedef int64_t qmckl_context_device;
