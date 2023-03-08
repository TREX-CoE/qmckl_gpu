// This files contains the C wrappers to be interfaced with Fortran

#include "qmckl_gpu.h"

//**********
// CONTEXT
//**********

qmckl_context_device
qmckl_context_touch_device_f(qmckl_context_device context) {
	return qmckl_context_touch_device(context);
}

qmckl_context_device qmckl_context_create_device_f(int *device_id) {
	return qmckl_context_create_device(*device_id);
}

qmckl_exit_code qmckl_context_destroy_device_f(qmckl_context_device context) {
	return qmckl_context_destroy_device(context);
}

//**********
// MEMORY
//**********

/* Allocs & frees */

void *qmckl_malloc_host_f(qmckl_context_device context,
						  const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host_f(qmckl_context_device context,
								  void *const ptr);

void *qmckl_malloc_device_f(qmckl_context_device context, size_t *size) {
	return qmckl_malloc_device(context, *size);
}

qmckl_exit_code qmckl_free_device_f(qmckl_context_device context,
									void *const ptr) {
	return qmckl_free_device(context, ptr);
}

/* Memcpys */

qmckl_exit_code qmckl_memcpy_H2D_f(qmckl_context_device context,
								   void *const dest, void *const src,
								   size_t *size) {
	return qmckl_memcpy_H2D(context, dest, src, *size);
}

qmckl_exit_code qmckl_memcpy_D2H_f(qmckl_context_device context,
								   void *const dest, void *const src,
								   size_t *size) {
	return qmckl_memcpy_D2H(context, dest, src, *size);
}

qmckl_exit_code qmckl_memcpy_D2D_f(qmckl_context_device context,
								   void *const dest, void *const src,
								   size_t *size) {
	return qmckl_memcpy_D2D(context, dest, src, *size);
}

//**********
// TREXIO
//**********

qmckl_exit_code qmckl_trexio_read_device_f(qmckl_context_device context,
										   char *file_name, int64_t *size_max) {
	return qmckl_trexio_read_device(context, file_name, *size_max);
}

qmckl_exit_code qmckl_set_electron_num_device_f(qmckl_context_device context,
												int64_t *up_num,
												int64_t *down_num) {
	return qmckl_set_electron_num_device(context, *up_num, *down_num);
}

qmckl_exit_code qmckl_set_electron_coord_device_f(qmckl_context_device context,
												  char *transp,
												  int64_t *walk_num,
												  double *coord,
												  int64_t *size_max) {
	return qmckl_set_electron_coord_device(context, *transp, *walk_num, coord,
										   *size_max);
}

qmckl_exit_code qmckl_set_point_device_f(qmckl_context_device context,
										 char *transp, int64_t *num,
										 double *coord, int64_t *size_max) {
	return qmckl_set_point_device(context, *transp, *num, coord, *size_max);
}

qmckl_exit_code qmckl_get_nucleus_num_device_f(qmckl_context_device context,
											   int64_t *num) {
	return qmckl_get_nucleus_num_device(context, num);
}

qmckl_exit_code qmckl_set_nucleus_num_device_f(qmckl_context_device context,
											   int64_t *num) {
	return qmckl_set_nucleus_num_device(context, *num);
}

qmckl_exit_code qmckl_set_nucleus_charge_device_f(qmckl_context_device context,
												  double *charge,
												  int64_t *size_max) {
	return qmckl_set_nucleus_charge_device(contextn charge, *size_max);
}

qmckl_exit_code qmckl_set_nucleus_coord_device_f(qmckl_context_device context,
												 char *transp, double *coord,
												 int64_t *size_max) {
	return qmckl_set_nucleus_coord_device(context, *transp, coord, *size_max);
}

qmckl_exit_code qmckl_get_ao_basis_ao_num_device_f(qmckl_context_device context,
												   int64_t *ao_num) {
	return qmckl_get_ao_basis_ao_num_device(context, ao_num);
}

qmckl_exit_code qmckl_set_ao_basis_type_device_f(qmckl_context_device context,
												 char *basis_type) {
	return qmckl_set_ao_basis_type_device(context, *basis_type);
}

qmckl_exit_code
qmckl_set_ao_basis_shell_num_device_f(qmckl_context_device context,
									  int64_t *shell_num) {
	return qmckl_set_ao_basis_shell_num_device(context, *shell_num);
}

qmckl_exit_code
qmckl_set_ao_basis_prim_num_device_f(qmckl_context_device context,
									 int64_t *prim_num) {
	return qmckl_set_ao_basis_prim_num_device(context, *prim_num);
}

qmckl_exit_code qmckl_set_ao_basis_coefficient_device_f(
	qmckl_context_device context, double *coefficient, int64_t *size_max) {
	return qmckl_set_ao_basis_coefficient_device(context, coefficient,
												 *size_max);
}

qmckl_exit_code qmckl_set_ao_basis_ao_num_device_f(qmckl_context_device context,
												   int64_t *ao_num) {
	return qmckl_set_ao_basis_ao_num_device(context, *ao_num);
}

qmckl_exit_code qmckl_set_ao_basis_nucleus_index_device_f(
	qmckl_context_device context, int64_t *nucleus_index, int64_t *size_max) {
	return qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index,
												   *size_max);
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num_device_f(qmckl_context_device context,
											  int64_t *nucleus_shell_num,
											  int64_t *size_max) {
	return qmckl_set_ao_basis_nucleus_shell_num_device(
		context, nucleus_shell_num, *size_max);
}

qmckl_exit_code qmckl_set_ao_basis_shell_ang_mom_device_f(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t *size_max) {
	return qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom,
												   *size_max);
}

qmckl_exit_code qmckl_set_ao_basis_shell_prim_num_device_f(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t *size_max) {
	return qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num,
													*size_max);
}

qmckl_exit_code qmckl_set_ao_basis_shell_prim_index_device_f(
	qmckl_context context, int64_t *shell_prim_index, int64_t *size_max) {
	qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index,
											   *size_max);
}

qmckl_exit_code qmckl_set_ao_basis_shell_factor_device_f(
	qmckl_context_device context, double *shell_factor, int64_t *size_max) {
	return qmckl_set_ao_basis_shell_factor_device(context, shell_factor,
												  *size_max);
}

qmckl_exit_code
qmckl_set_ao_basis_exponent_device_f(qmckl_context_device context,
									 double *exponent, int64_t *size_max) {
	return qmckl_set_ao_basis_exponent_device(context, exponent, *size_max);
}

qmckl_exit_code qmckl_set_ao_basis_prim_factor_device_f(qmckl_context context,
														double *prim_factor,
														int64_t *size_max) {
	return qmckl_set_ao_basis_prim_factor_device(context, prim_factor,
												 *size_max);
}

qmckl_exit_code
qmckl_set_ao_basis_ao_factor_device_f(qmckl_context_device context,
									  double *ao_factor, int64_t *size_max) {
	return qmckl_set_ao_basis_ao_factor_device(context, ao_factor, *size_max);
}

qmckl_exit_code
qmckl_get_ao_basis_shell_num_device_f(qmckl_context_device context,
									  int64_t *ao_num) {
	return qmckl_get_ao_basis_shell_num_device(context, ao_num);
}

qmckl_exit_code
qmckl_get_ao_basis_prim_num_device_f(qmckl_context_device context,
									 int64_t *prim_num) {
	return qmckl_get_ao_basis_prim_num_device(context, prim_num);
}

qmckl_exit_code qmckl_get_ao_basis_type_device_f(qmckl_context_device context,
												 char *type) {
	return qmckl_get_ao_basis_type_device(context, type);
}

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num_device_f(qmckl_context_device context,
											  int64_t *nucleus_shell_num,
											  int64_t *nucl_num) {
	return qmckl_get_ao_basis_nucleus_shell_num_device(
		context, nucleus_shell_num, *nucl_num);
}

qmckl_exit_code qmckl_get_ao_basis_shell_ang_mom_device_f(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t *shell_num) {
	return qmckl_get_ao_basis_shell_ang_mom_device(context, shell_ang_mom,
												   *shell_num);
}

qmckl_exit_code qmckl_get_ao_basis_shell_factor_device_f(
	qmckl_context_device context, double *shell_factor, int64_t *shell_num) {
	return qmckl_get_ao_basis_shell_factor_device(context, shell_factor,
												  *shell_num);
}

qmckl_exit_code qmckl_get_ao_basis_shell_prim_num_device_f(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t *shell_num) {
	return qmckl_get_ao_basis_shell_prim_num_device(context, shell_prim_num,
													*shell_num);
}

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index_device_f(qmckl_context_device context,
											 int64_t *nucleus_shell_num,
											 int64_t *shell_num) {
	return qmckl_get_ao_basis_shell_prim_index_device(
		context, nucleus_shell_num, *shell_num);
}

qmckl_exit_code
qmckl_get_ao_basis_exponent_device_f(qmckl_context_device context,
									 double *exponent, int64_t *prim_num) {
	return qmckl_get_ao_basis_exponent_device(context, exponent, *prim_num);
}

qmckl_exit_code qmckl_get_ao_basis_coefficient_device_f(
	qmckl_context_device context, double *coefficient, int64_t *prim_num) {
	return qmckl_get_ao_basis_coefficient_device(context, coefficient,
												 *prim_num);
}

qmckl_exit_code qmckl_get_ao_basis_prim_factor_device_f(
	qmckl_context_device context, double *prim_factor, int64_t *prim_num) {
	return qmckl_get_ao_basis_prim_factor_device(context, prim_factor,
												 *prim_num);
}

qmckl_exit_code
qmckl_get_ao_basis_ao_factor_device_f(qmckl_context_device context,
									  double *ao_factor, int64_t *ao_num) {
	return qmckl_get_ao_basis_ao_factor_device(context, ao_factor, *ao_num);
}

qmckl_exit_code qmckl_get_ao_basis_nucleus_index_device_f(
	qmckl_context_device context, int64_t *nucleus_index, int64_t *nucl_num) {
	return qmckl_get_ao_basis_nucleus_index_device(context, nucleus_index,
												   *nucl_num);
}

qmckl_exit_code qmckl_set_mo_basis_coefficient_device_f(qmckl_context context,
														double *coefficient) {
	return qmckl_set_mo_basis_coefficient_device(context, coefficient);
}

//**********
// AO
//**********

qmckl_exit_code qmckl_get_ao_basis_ao_num_device_f(qmckl_context_device context,
												   int64_t *ao_num) {
	return qmckl_get_ao_basis_ao_num_device(context, ao_num);
}

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device_f(qmckl_context_device context,
												   double *ao_vgl,
												   int64_t *size_max) {
	return qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl, *size_max);
}

qmckl_exit_code
qmckl_get_ao_basis_ao_value_device_f(qmckl_context_device context,
									 double *ao_vgl, int64_t *size_max) {
	return qmckl_get_ao_basis_ao_value_device(context, ao_vgl, *size_max);
}

//**********
// MO
//**********

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_device_f(qmckl_context context,
												   double *const mo_vgl,
												   const int64_t *size_max) {
	return qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl, *size_max);
}

qmckl_exit_code
qmckl_get_mo_basis_mo_value_device_f(qmckl_context_device context,
									 double *mo_value, int64_t *size_max) {
	return qmckl_get_mo_basis_mo_value_device(context, mo_vgl, *size_max);
}

qmckl_exit_code qmckl_get_mo_basis_mo_value_inplace_device_f(
	qmckl_context_device context, double *mo_value, int64_t *size_max) {
	return qmckl_get_mo_basis_mo_value_inplace_device(context, mo_value,
													  *size_max);
}

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_inplace_device_f(
	qmckl_context context, double *const mo_vgl, const int64_t *size_max) {
	return qmckl_get_mo_basis_mo_vgl_inplace_device(context, mo_vgl, *size_max);
}
