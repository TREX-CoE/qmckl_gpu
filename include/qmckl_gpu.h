// This is the header file meant to be included by the users.
// It contains prototypes for all GPU public functions, and definition of
// the _device context variant.

#pragma once

#include <qmckl.h>

typedef int64_t qmckl_context_device;
typedef int64_t qmckl_memory_info_struct;

//**********
// CONTEXT
//**********

qmckl_context_device qmckl_context_touch_device(qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);

qmckl_exit_code qmckl_context_destroy_device(qmckl_context_device context);

//**********
// ELECTRON
//**********

/* Allocs & frees */

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_device(qmckl_context_device context, size_t size);

qmckl_exit_code qmckl_free_device(qmckl_context_device context,
								  void *const ptr);

/* Memcpys */

qmckl_exit_code qmckl_memcpy_H2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code qmckl_memcpy_D2H(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code qmckl_memcpy_D2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);

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

bool qmckl_mo_basis_select_mo_device(qmckl_context_device context,
									 int32_t *keep,
									 int64_t size_max);

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
