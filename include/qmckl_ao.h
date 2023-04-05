#pragma once

// This file contains prototypes for the AO functions

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_basic_functions.h"

#include "qmckl_context.h"
#include "qmckl_memory.h"

qmckl_exit_code_device qmckl_compute_ao_basis_shell_gaussian_vgl_device(
	qmckl_context_device context, int prim_num, int shell_num, int point_num,
	int nucl_num, int64_t *nucleus_shell_num, int64_t *nucleus_index,
	double *nucleus_range, int64_t *shell_prim_index, int64_t *shell_prim_num,
	double *coord, double *nucl_coord, double *expo, double *coef_normalized,
	double *shell_vgl);

qmckl_exit_code_device qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_vgl);

qmckl_exit_code_device qmckl_compute_ao_value_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_value);

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_value_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
								 double *const ao_vgl, const int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_value_device(
	qmckl_context_device context, double *const ao_vgl, const int64_t size_max);
qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(const qmckl_context_device context,
									int64_t *const shell_num);
qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(const qmckl_context_device context,
								   int64_t *const prim_num);
qmckl_exit_code_device
qmckl_get_ao_basis_type_device(const qmckl_context_device context,
							   char *const basis_type);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);
qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_finalize_ao_basis_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_set_ao_basis_type_device(qmckl_context_device context, char basis_type);
qmckl_exit_code_device
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num);
qmckl_exit_code_device
qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);
qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_set_ao_basis_ao_num_device(qmckl_context_device context, int64_t ao_num);
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
qmckl_exit_code_device
qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t *ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t *prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_type_device(qmckl_context_device context, char *type);

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t nucl_num);
qmckl_exit_code_device qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t nucl_num);

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

qmckl_exit_code_device
qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t ao_num);

bool qmckl_ao_basis_provided(qmckl_context_device context);
