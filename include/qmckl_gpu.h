// This is the header file meant to be included by the users.
// It contains prototypes for all GPU public functions, and definition of
// the _device context variant.

#pragma once

#include <qmckl.h>

typedef int64_t qmckl_context_device;

//**********
// CONTEXT
//**********

qmckl_context_device
qmckl_context_touch_device(const qmckl_context_device context);
qmckl_context_device
qmckl_context_touch_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);
qmckl_context_device qmckl_context_create_device();

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context);

//**********
// TREXIO
//**********

qmckl_exit_code qmckl_trexio_read_device(const qmckl_context_device context,
											 const char *file_name,
											 const int64_t size_max);

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

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(
	qmckl_context_device context, double *const ao_vgl, const int64_t size_max);

//**********
// JASTROW
//**********

// TODO Jastrow has been reworked recently

//**********
// ELECTRON
//**********

qmckl_exit_code qmckl_set_electron_coord_device(
	qmckl_context_device context, const char transp, const int64_t walk_num,
	const double *coord, const int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
										   const char transp, const int64_t num,
										   const double *coord,
										   const int64_t size_max);
