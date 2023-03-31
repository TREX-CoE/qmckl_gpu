// This is the header file meant to be included by the users.
// It contains prototypes for all GPU public functions, and definition of
// the _device context variant.

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifdef HAVE_CUBLAS 
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

#ifdef HAVE_CUSPARSE 
#include <cuda_runtime_api.h>
#include <cusparse_v2.h>
#endif

/* CPU forward decls */
typedef int32_t qmckl_exit_code_device;
typedef int64_t qmckl_context_device;
typedef struct qmckl_memory_info_struct_device qmckl_memory_info_struct_device;
typedef struct qmckl_memory_struct_device qmckl_memory_struct_device;

/* QMCKL DEFINE */

#define VALID_TAG_DEVICE 0xBEEFFACE
#define INVALID_TAG_DEVICE 0xDEADBEEF

#define QMCKL_DEFAULT_PRECISION_DEVICE 53
#define QMCKL_DEFAULT_RANGE_DEVICE 11

#define QMCKL_MAX_FUN_LEN_DEVICE 256
#define QMCKL_MAX_MSG_LEN_DEVICE 1024

#define qmckl_mat_device(m, i, j) m.data[(i) + (j)*m.size[0]]

#define QMCKL_SUCCESS_DEVICE ((qmckl_exit_code_device)0)
#define QMCKL_INVALID_ARG_1_DEVICE ((qmckl_exit_code_device)1)
#define QMCKL_INVALID_ARG_2_DEVICE ((qmckl_exit_code_device)2)
#define QMCKL_INVALID_ARG_3_DEVICE ((qmckl_exit_code_device)3)
#define QMCKL_INVALID_ARG_4_DEVICE ((qmckl_exit_code_device)4)
#define QMCKL_INVALID_ARG_5_DEVICE ((qmckl_exit_code_device)5)
#define QMCKL_INVALID_ARG_6_DEVICE ((qmckl_exit_code_device)6)
#define QMCKL_INVALID_ARG_7_DEVICE ((qmckl_exit_code_device)7)
#define QMCKL_INVALID_ARG_8_DEVICE ((qmckl_exit_code_device)8)
#define QMCKL_INVALID_ARG_9_DEVICE ((qmckl_exit_code_device)9)
#define QMCKL_INVALID_ARG_10_DEVICE ((qmckl_exit_code_device)10)
#define QMCKL_INVALID_ARG_11_DEVICE ((qmckl_exit_code_device)11)
#define QMCKL_INVALID_ARG_12_DEVICE ((qmckl_exit_code_device)12)
#define QMCKL_INVALID_ARG_13_DEVICE ((qmckl_exit_code_device)13)
#define QMCKL_INVALID_ARG_14_DEVICE ((qmckl_exit_code_device)14)
#define QMCKL_INVALID_ARG_15_DEVICE ((qmckl_exit_code_device)15)
#define QMCKL_INVALID_ARG_16_DEVICE ((qmckl_exit_code_device)16)
#define QMCKL_INVALID_ARG_17_DEVICE ((qmckl_exit_code_device)17)
#define QMCKL_INVALID_ARG_18_DEVICE ((qmckl_exit_code_device)18)
#define QMCKL_INVALID_ARG_19_DEVICE ((qmckl_exit_code_device)19)
#define QMCKL_INVALID_ARG_20_DEVICE ((qmckl_exit_code_device)20)
#define QMCKL_FAILURE_DEVICE ((qmckl_exit_code_device)101)
#define QMCKL_ERRNO_DEVICE ((qmckl_exit_code_device)102)
#define QMCKL_INVALID_CONTEXT_DEVICE ((qmckl_exit_code_device)103)
#define QMCKL_ALLOCATION_FAILED_DEVICE ((qmckl_exit_code_device)104)
#define QMCKL_DEALLOCATION_FAILED_DEVICE ((qmckl_exit_code_device)105)
#define QMCKL_NOT_PROVIDED_DEVICE ((qmckl_exit_code_device)106)
#define QMCKL_OUT_OF_BOUNDS_DEVICE ((qmckl_exit_code_device)107)
#define QMCKL_ALREADY_SET_DEVICE ((qmckl_exit_code_device)108)
#define QMCKL_INVALID_EXIT_CODE_DEVICE ((qmckl_exit_code_device)109)

//**********
// CONTEXT
//**********

qmckl_exit_code_device qmckl_context_touch_device(qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);

qmckl_exit_code_device qmckl_context_destroy_device(qmckl_context_device context);

//**********
// MEMORY
//**********

typedef struct qmckl_memory_info_struct qmckl_memory_info_struct;

/* Allocs & frees */
void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code_device qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_device(qmckl_context_device context, size_t size);

qmckl_exit_code_device qmckl_free_device(qmckl_context_device context, void *ptr);

/* Memcpys */

qmckl_exit_code_device qmckl_memcpy_H2D(qmckl_context_device context, void *dest,
								 void *src, size_t size);
qmckl_exit_code_device qmckl_memcpy_D2H(qmckl_context_device context, void *dest,
								 void *src, size_t size);
qmckl_exit_code_device qmckl_memcpy_D2D(qmckl_context_device context, void *dest,
								 void *src, size_t size);

// *********
// ERROR MANAGEMENT
// *********

const char *qmckl_string_of_error(const qmckl_exit_code_device error);

//**********
// TREXIO
//**********

qmckl_exit_code_device qmckl_trexio_read_device(qmckl_context_device context,
										 char *file_name, int64_t size_max);

qmckl_exit_code_device qmckl_set_electron_num_device(qmckl_context_device context,
											  int64_t up_num, int64_t down_num);

qmckl_exit_code_device qmckl_set_electron_coord_device(qmckl_context_device context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max);

qmckl_exit_code_device qmckl_set_point_device(qmckl_context_device context,
									   char transp, int64_t num, double *coord,
									   int64_t size_max);

qmckl_exit_code_device qmckl_get_nucleus_num_device(qmckl_context_device context,
											 int64_t *num);

qmckl_exit_code_device qmckl_set_nucleus_num_device(qmckl_context_device context,
											 int64_t num);

qmckl_exit_code_device qmckl_set_nucleus_charge_device(qmckl_context_device context,
												double *charge,
												int64_t size_max);
qmckl_exit_code_device qmckl_set_nucleus_coord_device(qmckl_context_device context,
											   char transp, double *coord,
											   int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code_device qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   char basis_type);
qmckl_exit_code_device
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num);
qmckl_exit_code_device qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t prim_num);

qmckl_exit_code_device qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code_device
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max);

qmckl_exit_code_device qmckl_set_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t ao_num);
qmckl_exit_code_device qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max);

qmckl_exit_code_device qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t size_max);

qmckl_exit_code_device qmckl_set_ao_basis_prim_factor_device(qmckl_context_device context,
													  double *prim_factor,
													  int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t *ao_num);

qmckl_exit_code_device qmckl_get_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t *prim_num);

qmckl_exit_code_device qmckl_get_ao_basis_type_device(qmckl_context_device context,
											   char *type);

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t nucl_num);

qmckl_exit_code_device qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t shell_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t shell_num);

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t shell_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_prim_index_device(qmckl_context_device context,
										   int64_t *nucleus_shell_num,
										   int64_t shell_num);

qmckl_exit_code_device qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t ao_num);

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t nucl_num);

qmckl_exit_code_device qmckl_set_mo_basis_coefficient_device(qmckl_context_device context,
													  double *coefficient);
qmckl_exit_code_device qmckl_set_mo_basis_mo_num_device(qmckl_context_device context,
												 int64_t mo_num);
bool qmckl_nucleus_provided(qmckl_context_device context);

//**********
// AO
//**********

qmckl_exit_code_device qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code_device qmckl_get_ao_basis_ao_vgl_inplace_offload(qmckl_context_device context,
														  double *ao_vgl,
														  int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *ao_vgl,
												 int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_value_device(qmckl_context_device context,
												   double *ao_vgl,
												   int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_vgl_acc_offload(qmckl_context_device context,
													  double *const ao_vgl,
													  const int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_vgl_inplace_acc_offload(
	qmckl_context_device context, double *const ao_vgl, const int64_t size_max);

bool qmckl_ao_basis_provided(qmckl_context_device context);

//**********
// MO
//**********

bool qmckl_mo_basis_select_mo_device(qmckl_context_device context,
									 int32_t *keep, int64_t size_max);

qmckl_exit_code_device qmckl_get_mo_basis_mo_num(const qmckl_context_device context,
										  int64_t *mo_num);

qmckl_exit_code_device qmckl_get_mo_basis_mo_vgl_device(qmckl_context_device context,
												 double *const mo_vgl,
												 const int64_t size_max);

qmckl_exit_code_device qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
												   double *mo_value,
												   int64_t size_max);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *mo_value, int64_t size_max);

qmckl_exit_code_device qmckl_get_mo_basis_mo_vgl_inplace_device(
	qmckl_context_device context, double *const mo_vgl, const int64_t size_max);

qmckl_exit_code_device qmckl_provide_mo_basis_mo_value_device(qmckl_context_device context);

qmckl_exit_code_device qmckl_provide_mo_basis_mo_vgl_device(qmckl_context_device context);

qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num, int64_t point_num,
	double *restrict coefficient_t, double *restrict ao_value,
	double *restrict mo_value);

qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num, int64_t point_num,
	double *restrict coefficient_t, double *restrict ao_vgl,
	double *restrict mo_vgl);

qmckl_exit_code_device qmckl_get_mo_basis_mo_vgl_acc_offload(qmckl_context_device context,
													  double *const mo_vgl,
													  const int64_t size_max);

qmckl_exit_code_device qmckl_get_mo_basis_mo_vgl_acc_offload_inplace(
	qmckl_context_device context, double *const mo_vgl, const int64_t size_max);

qmckl_exit_code_device qmckl_get_mo_basis_mo_num_device(const qmckl_context_device context,
												 int64_t *mo_num);

bool qmckl_mo_basis_provided(qmckl_context_device context);

//**********
// ELECTRON
//**********

qmckl_exit_code_device qmckl_set_electron_coord_device(qmckl_context_device context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max);

qmckl_exit_code_device qmckl_set_point_device(qmckl_context_device context,
									   char transp, int64_t num, double *coord,
									   int64_t size_max);

//****************************
// SHERMAN-MORRISON & WOODBURY
//****************************

// TODO In the future, we could to generate this header on the fly so it
// contains exactly the functions that are enabled (no need for the end user to
// worry about preprocessor)

#ifdef HAVE_CUBLAS
qmckl_exit_code_device
qmckl_woodbury_kxk(const qmckl_context_device context,
                   cublasHandle_t b_handle, cusolverDnHandle_t s_handle,
				   const uint64_t Lds, const uint64_t Dim,
				   const uint64_t N_updates, const double *__restrict Updates,
				   const uint64_t *__restrict Updates_index,
				   const double breakdown, double *__restrict Slater_inv,
				   double *__restrict determinant);
#endif
