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

qmckl_context_device
qmckl_context_touch_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context);

//**********
// ELECTRON
//**********

/* Allocs & frees */

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_device(qmckl_context_device context,
						  const qmckl_memory_info_struct info);

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

qmckl_exit_code qmckl_trexio_read_device(const qmckl_context_device context,
										 const char *file_name,
										 const int64_t size_max);

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  const int64_t up_num,
											  const int64_t down_num);

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  const int64_t up_num,
											  const int64_t down_num);

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context_device context,
												const char transp,
												const int64_t walk_num,
												const double *coord,
												const int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   const char transp, const int64_t num,
									   const double *coord,
									   const int64_t size_max);

qmckl_exit_code qmckl_get_nucleus_num_device(const qmckl_context_device context,
											 int64_t *const num);
qmckl_exit_code qmckl_get_nucleus_num_device(const qmckl_context_device context,
											 int64_t *const num);

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 const int64_t num);

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 const int64_t num);

qmckl_exit_code qmckl_set_nucleus_charge_device(qmckl_context_device context,
												const double *charge,
												const int64_t size_max);
qmckl_exit_code qmckl_set_nucleus_coord_device(qmckl_context_device context,
											   const char transp,
											   const double *coord,
											   const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);

qmckl_exit_code qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   const char basis_type);
qmckl_exit_code
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									const int64_t shell_num);
qmckl_exit_code qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   const int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);

qmckl_exit_code qmckl_set_ao_basis_ao_num_device(qmckl_context_device context,
												 const int64_t ao_num);
qmckl_exit_code
qmckl_set_ao_basis_nucleus_index_device(qmckl_context_device context,
										const int64_t *nucleus_index,
										const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											const int64_t *nucleus_shell_num,
											const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_shell_ang_mom_device(qmckl_context_device context,
										const int32_t *shell_ang_mom,
										const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_shell_prim_num_device(qmckl_context_device context,
										 const int64_t *shell_prim_num,
										 const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_shell_prim_index_device(qmckl_context context,
										   const int64_t *shell_prim_index,
										   const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   const double *shell_factor,
									   const int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
												   const double *exponent,
												   const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  const double *coefficient,
									  const int64_t size_max);
qmckl_exit_code qmckl_set_ao_basis_prim_factor_device(qmckl_context context,
													  const double *prim_factor,
													  const int64_t size_max);
qmckl_exit_code
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									const double *ao_factor,
									const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_num_device(const qmckl_context_device context,
									int64_t *const ao_num);

qmckl_exit_code
qmckl_get_ao_basis_prim_num_device(const qmckl_context_device context,
								   int64_t *const prim_num);

qmckl_exit_code
qmckl_get_ao_basis_type_device(const qmckl_context_device context,
							   char *const type);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num_device(const qmckl_context_device context,
											int64_t *const nucleus_shell_num,
											int64_t nucl_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_ang_mom_device(const qmckl_context_device context,
										double *const shell_ang_mom,
										int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_factor_device(const qmckl_context_device context,
									   int32_t *const shell_factor,
									   int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_num_device(const qmckl_context_device context,
										 double *const shell_prim_num,
										 int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index_device(const qmckl_context_device context,
										   int64_t *const nucleus_shell_num,
										   int64_t shell_num);

qmckl_exit_code
qmckl_get_ao_basis_exponent_device(const qmckl_context_device context,
								   double *const exponent, int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_coefficient_device(const qmckl_context_device context,
									  double *const coefficient,
									  int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_prim_factor_device(const qmckl_context_device context,
									  double *const prim_factor,
									  int64_t prim_num);

qmckl_exit_code
qmckl_get_ao_basis_ao_factor_device(const qmckl_context_device context,
									double *const ao_factor, int64_t ao_num);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_index_device(const qmckl_context_device context,
										int64_t *const nucleus_index,
										int64_t nucl_num);


//**********
// AO
//**********

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_offload(
	qmckl_context context, double *const ao_vgl, const int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *const ao_vgl,
												 const int64_t size_max);

//**********
// JASTROW
//**********

// TODO Jastrow has been reworked recently

//**********
// ELECTRON
//**********

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context_device context,
												const char transp,
												const int64_t walk_num,
												const double *coord,
												const int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   const char transp, const int64_t num,
									   const double *coord,
									   const int64_t size_max);
