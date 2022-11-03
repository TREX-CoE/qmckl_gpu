// This is the header file meant to be included by the users.
// It contains prototypes for all GPU functions, and definition of
// _device variants fo the structure datatypes.

#pragma once

#include <qmckl.h>

typedef int64_t qmckl_context_device;

//**********
// CONTEXT
//**********

qmckl_context_device
qmckl_context_check_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device();

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context);

//**********
// TREXIO
//**********

qmckl_exit_code qmckl_trexio_read_device(const qmckl_context_device context,
                                         const char *file_name,
                                         const int64_t size_max, int device_id);

//**********
// AO
//**********

// qmckl_ao_openmp.c
qmckl_exit_code qmckl_get_ao_basis_ao_vgl_omp_offload(qmckl_context context,
                                                      double *const ao_vgl,
                                                      const int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_omp_offload(
    qmckl_context context, double *const ao_vgl, const int64_t size_max);

// qmckl_ao_acc.c
qmckl_exit_code qmckl_get_ao_basis_ao_vgl_acc_offload(qmckl_context context,
                                                      double *const ao_vgl,
                                                      const int64_t size_max);

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_acc_offload(
    qmckl_context context, double *const ao_vgl, const int64_t size_max);

// qmckl_ao_device.c
qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
                                                 double *const ao_vgl,
                                                 const int64_t size_max,
                                                 int device_id);

//**********
// JASTROW
//**********

// qmckl_jastrow_openmp.c
qmckl_exit_code qmckl_get_jastrow_tmp_c_omp_offload(qmckl_context context,
                                                    double *const tmp_c);

qmckl_exit_code qmckl_get_jastrow_dtmp_c_omp_offload(qmckl_context context,
                                                     double *const dtmp_c);

//**********
// ELECTRON
//**********

qmckl_exit_code
qmckl_set_electron_coord_device(qmckl_context_device context,
                                const char transp,
                                const int64_t walk_num,
                                const double* coord,
                                const int64_t size_max,
                                int device_id);
