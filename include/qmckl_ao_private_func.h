#ifndef QMCKL_AO_HPF
#define QMCKL_AO_HPF
#include "qmckl_blas_private_type.h"



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_ao_basis(qmckl_context context);

/* After initialization */

/* When the basis set is completely entered, extra data structures may be */
/* computed to accelerate the calculations. The primitives within each */
/* contraction are sorted in ascending order of their exponents, such */
/* that as soon as a primitive is zero all the following functions */
/* vanish. Also, it is possible to compute a nuclear radius beyond which */
/* all the primitives are zero up to the numerical accuracy defined in */
/* the context. */


qmckl_exit_code qmckl_finalize_basis (qmckl_context context);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_finalize_basis_device (qmckl_context context, int device_id);
#endif

#ifdef HAVE_HPC
qmckl_exit_code qmckl_finalize_basis_hpc (qmckl_context context);
#endif

#if defined(HAVE_HPC) && defined(HAVE_DEVICE_POINTERS)
qmckl_exit_code qmckl_finalize_basis_hpc_device (qmckl_context context, int device_pointer);
#endif

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ao_basis_primitive_vgl(qmckl_context context);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ao_basis_shell_vgl(qmckl_context context);

/* Interfaces */
/*  #  #+CALL: generate_c_header(table=qmckl_ao_value_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_value")) */
/*  #  (Commented because the header needs to go into h_private_func) */


qmckl_exit_code qmckl_compute_ao_value_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const double* coord,
      const double* nucl_coord,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const double* nucleus_range,
      const int32_t* nucleus_max_ang_mom,
      const int32_t* shell_ang_mom,
      const double* ao_factor,
      const double* shell_vgl,
      double* const ao_value );

#ifdef HAVE_HPC
    qmckl_exit_code qmckl_compute_ao_value_hpc_gaussian (
          const qmckl_context context,
          const int64_t ao_num,
          const int64_t shell_num,
          const int32_t* prim_num_per_nucleus,
          const int64_t point_num,
          const int64_t nucl_num,
          const double* coord,
          const double* nucl_coord,
          const int64_t* nucleus_index,
          const int64_t* nucleus_shell_num,
          const double* nucleus_range,
          const int32_t* nucleus_max_ang_mom,
          const int32_t* shell_ang_mom,
          const double* ao_factor,
          const qmckl_matrix expo_per_nucleus,
          const qmckl_tensor coef_per_nucleus,
          double* const ao_value );
#endif

/* Provide                                                       :noexport: */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_value(qmckl_context context);

/* Device_pointers */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_value_device(qmckl_context context, int64_t device_id);

qmckl_exit_code qmckl_provide_ao_value(qmckl_context context);

/* Device pointers version */



#ifdef HAVE_DEVICE_POINTERS

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
#endif

/* Interfaces */
/*  #  #+CALL: generate_c_header(table=qmckl_ao_vgl_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_vgl")) */
/*  #  (Commented because the header needs to go into h_private_func) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_vgl_doc (
                                          const qmckl_context context,
                                          const int64_t ao_num,
                                          const int64_t shell_num,
                                          const int64_t point_num,
                                          const int64_t nucl_num,
                                          const double* coord,
                                          const double* nucl_coord,
                                          const int64_t* nucleus_index,
                                          const int64_t* nucleus_shell_num,
                                          const double* nucleus_range,
                                          const int32_t* nucleus_max_ang_mom,
                                          const int32_t* shell_ang_mom,
                                          const double* ao_factor,
                                          const double* shell_vgl,
                                          double* const ao_vgl );

#ifdef HAVE_OPENMP_OFFLOAD
qmckl_exit_code qmckl_compute_ao_vgl_gaussian_omp_offload (
                                                   const qmckl_context context,
                                                   const int64_t ao_num,
                                                   const int64_t shell_num,
                                                   const int32_t* prim_num_per_nucleus,
                                                   const int64_t point_num,
                                                   const int64_t nucl_num,
                                                   const double* coord,
                                                   const double* nucl_coord,
                                                   const int64_t* nucleus_index,
                                                   const int64_t* nucleus_shell_num,
                                                   const double* nucleus_range,
                                                   const int32_t* nucleus_max_ang_mom,
                                                   const int32_t* shell_ang_mom,
                                                   const double* ao_factor,
                                                   const qmckl_matrix expo_per_nucleus,
                                                   const qmckl_tensor coef_per_nucleus,
                                                   double* const ao_vgl );
#elif HAVE_HPC

qmckl_exit_code qmckl_compute_ao_vgl_hpc_gaussian (
                                                           const qmckl_context context,
                                                           const int64_t ao_num,
                                                           const int64_t shell_num,
                                                           const int32_t* prim_num_per_nucleus,
                                                           const int64_t point_num,
                                                           const int64_t nucl_num,
                                                           const double* coord,
                                                           const double* nucl_coord,
                                                           const int64_t* nucleus_index,
                                                           const int64_t* nucleus_shell_num,
                                                           const double* nucleus_range,
                                                           const int32_t* nucleus_max_ang_mom,
                                                           const int32_t* shell_ang_mom,
                                                           const double* ao_factor,
                                                           const qmckl_matrix expo_per_nucleus,
                                                           const qmckl_tensor coef_per_nucleus,
                                                           double* const ao_vgl );
#endif

/* Provide                                                       :noexport: */

/* #+CALL: write_provider_header( group="ao_basis", data="ao_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_vgl(qmckl_context context);

/* DEVICE POITERS */


qmckl_exit_code qmckl_provide_ao_basis_ao_vgl_device(qmckl_context context, int64_t const device_id);

qmckl_exit_code qmckl_provide_ao_vgl(qmckl_context context);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_ao_vgl_device(qmckl_context context, int device_id);
#endif

#endif
