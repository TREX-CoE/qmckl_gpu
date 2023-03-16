#pragma once

// This file contains prototypes for the MO functions

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_device(qmckl_context context,
												 double *const mo_vgl,
												 const int64_t size_max);

qmckl_exit_code qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
												   double *mo_value,
												   int64_t size_max);

qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *mo_value, int64_t size_max);

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
