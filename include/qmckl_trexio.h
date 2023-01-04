#pragma once

// This file contains functions prototypes for functions initializing the
// context from a TREXIO file

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <trexio.h>

#include <omp.h>

#include <qmckl.h>
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"

// Prototype for standard QMCkl function
trexio_t *qmckl_trexio_open_X(const char *file_name, qmckl_exit_code *rc);

// Prototypes for our new device variants

//**********
// ELECTRON/POINT GETTERS/SETTERS
//**********

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

//**********
// NUCLEUS GETTERS/SETTERS
//**********

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
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context);
qmckl_exit_code
qmckl_finalize_nucleus_basis_device(qmckl_context_device context);

//**********
// AO GETTERS/SETTERS
//**********

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);

qmckl_exit_code
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context);

qmckl_exit_code qmckl_finalize_ao_basis_device(qmckl_context_device context);

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

//**********
// MO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_finalize_mo_basis_device(qmckl_context_device context);
qmckl_exit_code qmckl_set_mo_basis_mo_num_device(qmckl_context_device context,
												 const int64_t mo_num);
qmckl_exit_code
qmckl_set_mo_basis_coefficient_device(qmckl_context context,
									  const double *coefficient);

//**********
// CONTEXT FILL
//**********

qmckl_exit_code qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
												   trexio_t *const file);
qmckl_exit_code qmckl_trexio_read_ao_X_device(qmckl_context context,
											  trexio_t *const file);
qmckl_exit_code qmckl_trexio_read_mo_X_device(qmckl_context_device context,
											  trexio_t *const file);
qmckl_exit_code qmckl_trexio_read_device(const qmckl_context_device context,
										 const char *file_name,
										 const int64_t size_max);
