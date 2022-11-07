#ifndef QMCKL_MO_HPF
#define QMCKL_MO_HPF
#include "qmckl_blas_private_type.h"

/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */

qmckl_exit_code qmckl_init_mo_basis(qmckl_context context);

/* When the basis set is completely entered, other data structures are */
/* computed to accelerate the calculations. */

qmckl_exit_code qmckl_finalize_mo_basis(qmckl_context context);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_finalize_mo_basis_device(qmckl_context context,
                                               int64_t device_id);
#endif

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_value" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value(qmckl_context context);

/* pointer device */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value_device(qmckl_context context,
                                                       const int64_t device_id);

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl(qmckl_context context);

/* Device pointers */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl_device(qmckl_context context,
                                                     int64_t device_id);

/* Mixte AO/MOxAO */

qmckl_exit_code
qmckl_provide_ao_mo_basis_ao_mo_vgl_device(qmckl_context context,
                                           int64_t device_id);

/* MOxAO in the same kernel */

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_ao_mo_basis_ao_mo_vgl_device(
    const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
    const int64_t shell_num, const int32_t *restrict prim_num_per_nucleus,
    const int64_t nucl_num, const double *restrict coord,
    const double *restrict nucl_coord, const int64_t *restrict nucleus_index,
    const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
    const int32_t *restrict nucleus_max_ang_mom,
    const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
    const qmckl_matrix expo_per_nucleus, const qmckl_tensor coef_per_nucleus,
    const int64_t point_num, const double *coefficient_t, double *ao_vgl,
    double *mo_vgl, const int64_t device_id);

#endif

#endif
