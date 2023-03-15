// This is the header file meant to be included by the users.
// It contains prototypes for all GPU public functions, and definition of
// the _device context variant.

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* CPU */
typedef int32_t qmckl_exit_code;
typedef int64_t qmckl_context ;
/* GPU */
typedef int64_t qmckl_context_device;
//typedef int64_t qmckl_memory_info_struct;


/* QMCKL DEFINE */

#define VALID_TAG   0xBEEFFACE
#define INVALID_TAG 0xDEADBEEF

#define  QMCKL_DEFAULT_PRECISION        53
#define  QMCKL_DEFAULT_RANGE            11

#define  QMCKL_MAX_FUN_LEN   256
#define  QMCKL_MAX_MSG_LEN  1024


#define qmckl_mat(m, i, j) m.data[(i) + (j)*m.size[0]]


#define  QMCKL_SUCCESS                  ((qmckl_exit_code) 0)
#define  QMCKL_INVALID_ARG_1            ((qmckl_exit_code) 1)
#define  QMCKL_INVALID_ARG_2            ((qmckl_exit_code) 2)
#define  QMCKL_INVALID_ARG_3            ((qmckl_exit_code) 3)
#define  QMCKL_INVALID_ARG_4            ((qmckl_exit_code) 4)
#define  QMCKL_INVALID_ARG_5            ((qmckl_exit_code) 5)
#define  QMCKL_INVALID_ARG_6            ((qmckl_exit_code) 6)
#define  QMCKL_INVALID_ARG_7            ((qmckl_exit_code) 7)
#define  QMCKL_INVALID_ARG_8            ((qmckl_exit_code) 8)
#define  QMCKL_INVALID_ARG_9            ((qmckl_exit_code) 9)
#define  QMCKL_INVALID_ARG_10           ((qmckl_exit_code) 10)
#define  QMCKL_INVALID_ARG_11           ((qmckl_exit_code) 11)
#define  QMCKL_INVALID_ARG_12           ((qmckl_exit_code) 12)
#define  QMCKL_INVALID_ARG_13           ((qmckl_exit_code) 13)
#define  QMCKL_INVALID_ARG_14           ((qmckl_exit_code) 14)
#define  QMCKL_INVALID_ARG_15           ((qmckl_exit_code) 15)
#define  QMCKL_INVALID_ARG_16           ((qmckl_exit_code) 16)
#define  QMCKL_INVALID_ARG_17           ((qmckl_exit_code) 17)
#define  QMCKL_INVALID_ARG_18           ((qmckl_exit_code) 18)
#define  QMCKL_INVALID_ARG_19           ((qmckl_exit_code) 19)
#define  QMCKL_INVALID_ARG_20           ((qmckl_exit_code) 20)
#define  QMCKL_FAILURE                  ((qmckl_exit_code) 101)
#define  QMCKL_ERRNO                    ((qmckl_exit_code) 102)
#define  QMCKL_INVALID_CONTEXT          ((qmckl_exit_code) 103)
#define  QMCKL_ALLOCATION_FAILED        ((qmckl_exit_code) 104)
#define  QMCKL_DEALLOCATION_FAILED      ((qmckl_exit_code) 105)
#define  QMCKL_NOT_PROVIDED             ((qmckl_exit_code) 106)
#define  QMCKL_OUT_OF_BOUNDS            ((qmckl_exit_code) 107)
#define  QMCKL_ALREADY_SET              ((qmckl_exit_code) 108)
#define  QMCKL_INVALID_EXIT_CODE        ((qmckl_exit_code) 109)





//**********
// CONTEXT
//**********

qmckl_exit_code qmckl_context_touch_device(qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);

qmckl_exit_code qmckl_context_destroy_device(qmckl_context_device context);

//**********
// ELECTRON
//**********

/* Allocs & frees */
/*
void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_device(qmckl_context_device context, size_t size);

qmckl_exit_code qmckl_free_device(qmckl_context_device context, void *ptr);
*/
/* Memcpys */
/*
qmckl_exit_code qmckl_memcpy_H2D(qmckl_context_device context, void *dest,
								 void *src, size_t size);
qmckl_exit_code qmckl_memcpy_D2H(qmckl_context_device context, void *dest,
								 void *src, size_t size);
qmckl_exit_code qmckl_memcpy_D2D(qmckl_context_device context, void *dest,
								 void *src, size_t size);
*/

// *********
// ERROR MANAGEMENT
// *********

const char*
qmckl_string_of_error (const qmckl_exit_code error);




//**********
// TREXIO
//**********

qmckl_exit_code qmckl_trexio_read_device(qmckl_context_device context,
										 char *file_name, int64_t size_max);

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  int64_t up_num, int64_t down_num);

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context_device context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   char transp, int64_t num, double *coord,
									   int64_t size_max);

qmckl_exit_code qmckl_get_nucleus_num_device(qmckl_context_device context,
											 int64_t *num);

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 int64_t num);

qmckl_exit_code qmckl_set_nucleus_charge_device(qmckl_context_device context,
												double *charge,
												int64_t size_max);
qmckl_exit_code qmckl_set_nucleus_coord_device(qmckl_context_device context,
											   char transp, double *coord,
											   int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   char basis_type);
qmckl_exit_code
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num);
qmckl_exit_code qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t prim_num);

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max);

qmckl_exit_code qmckl_set_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t ao_num);
qmckl_exit_code qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context context, int64_t *shell_prim_index, int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max);

qmckl_exit_code qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t size_max);

qmckl_exit_code qmckl_set_ao_basis_prim_factor_device(qmckl_context context,
													  double *prim_factor,
													  int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t *ao_num);

qmckl_exit_code qmckl_get_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t *prim_num);

qmckl_exit_code qmckl_get_ao_basis_type_device(qmckl_context_device context,
											   char *type);

qmckl_exit_code qmckl_get_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t nucl_num);

qmckl_exit_code qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t shell_num);

qmckl_exit_code qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index_device(qmckl_context_device context,
										   int64_t *nucleus_shell_num,
										   int64_t shell_num);

qmckl_exit_code qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t ao_num);

qmckl_exit_code qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t nucl_num);

qmckl_exit_code qmckl_set_mo_basis_coefficient_device(qmckl_context context,
													  double *coefficient);

//**********
// AO
//**********

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_offload(qmckl_context context,
														  double *ao_vgl,
														  int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *ao_vgl,
												 int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_value_device(qmckl_context_device context,
												   double *ao_vgl,
												   int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_acc_offload(qmckl_context context,
													  double *const ao_vgl,
													  const int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_acc_offload(
	qmckl_context context, double *const ao_vgl, const int64_t size_max);

//**********
// MO
//**********

qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num);



qmckl_exit_code qmckl_get_mo_basis_mo_vgl_device(qmckl_context context,
												 double *const mo_vgl,
												 const int64_t size_max);

qmckl_exit_code qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
												   double *mo_value,
												   int64_t size_max);

qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *mo_value, int64_t size_max);

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_inplace_device(
	qmckl_context context, double *const mo_vgl, const int64_t size_max);

qmckl_exit_code qmckl_provide_mo_basis_mo_value_device(qmckl_context context);

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl_device(qmckl_context context);

qmckl_exit_code qmckl_compute_mo_basis_mo_value_device(
	qmckl_context context, int64_t ao_num, int64_t mo_num, int64_t point_num,
	double *restrict coefficient_t, double *restrict ao_value,
	double *restrict mo_value);

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_device(
	qmckl_context context, int64_t ao_num, int64_t mo_num, int64_t point_num,
	double *restrict coefficient_t, double *restrict ao_vgl,
	double *restrict mo_vgl);

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_acc_offload(qmckl_context context,
													  double *const mo_vgl,
													  const int64_t size_max);

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_acc_offload_inplace(
	qmckl_context context, double *const mo_vgl, const int64_t size_max);

//**********
// ELECTRON
//**********

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context_device context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   char transp, int64_t num, double *coord,
									   int64_t size_max);


//****************************
// SHERMAN-MORRISON & WOODBURY
//****************************

qmckl_exit_code qmckl_woodbury_kxk(
    cublasHandle_t b_handle,
    cusolverDnHandle_t s_handle,
    const uint64_t Lds,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* __restrict Updates,
    const uint64_t* __restrict Updates_index,
    const double breakdown,
    double* __restrict Slater_inv,
    double* __restrict determinant);

