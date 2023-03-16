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

#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"

// Prototype for standard QMCkl function
trexio_t *qmckl_trexio_open_X(char *file_name, qmckl_exit_code *rc);

// Prototypes for our new device variants

//**********
// ELECTRON/POINT GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  int64_t up_num, int64_t down_num);

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  int64_t up_num, int64_t down_num);

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context_device context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max);

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   char transp, int64_t num, double *coord,
									   int64_t size_max);

//**********
// NUCLEUS GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_get_nucleus_num_device(qmckl_context_device context,
											 int64_t *num);
qmckl_exit_code qmckl_get_nucleus_num_device(qmckl_context_device context,
											 int64_t *num);

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 int64_t num);

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 int64_t num);

qmckl_exit_code qmckl_set_nucleus_charge_device(qmckl_context_device context,
												double *charge,
												int64_t size_max);
qmckl_exit_code qmckl_set_nucleus_coord_device(qmckl_context_device context,
											   char transp, double *coord,
											   int64_t size_max);

qmckl_exit_code
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context);
qmckl_exit_code
qmckl_finalize_nucleus_basis_device(qmckl_context_device context);

//**********
// AO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);
qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

qmckl_exit_code
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context);

qmckl_exit_code qmckl_finalize_ao_basis_device(qmckl_context_device context);

qmckl_exit_code qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   char basis_type);
qmckl_exit_code
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num);
qmckl_exit_code qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t prim_num);

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);
qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num);

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
qmckl_exit_code
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max);
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
qmckl_exit_code qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t nucl_num);

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

//**********
// MO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_finalize_mo_basis_device(qmckl_context_device context);
qmckl_exit_code qmckl_set_mo_basis_mo_num_device(qmckl_context_device context,
												 int64_t mo_num);
qmckl_exit_code qmckl_set_mo_basis_coefficient_device(qmckl_context context,
													  double *coefficient);

//**********
// CONTEXT FILL
//**********

qmckl_exit_code qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
												   trexio_t *file);
qmckl_exit_code qmckl_trexio_read_ao_X_device(qmckl_context context,
											  trexio_t *file);
qmckl_exit_code qmckl_trexio_read_mo_X_device(qmckl_context_device context,
											  trexio_t *file);
qmckl_exit_code qmckl_trexio_read_device(qmckl_context_device context,
										 char *file_name, int64_t size_max);
