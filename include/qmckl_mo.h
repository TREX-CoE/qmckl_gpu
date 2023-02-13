#pragma once

// This file contains prototypes for the MO functions

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <qmckl.h>

#include "qmckl_ao_private_func.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_blas_private_type.h"
#include "qmckl_context_private_type.h"
#include "qmckl_error_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_memory_private_type.h"

#include "qmckl_context.h"
#include "qmckl_memory.h"

qmckl_exit_code qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
												   double *const mo_value,
												   const int64_t size_max);

qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *const mo_value,
										   const int64_t size_max);

qmckl_exit_code qmckl_provide_mo_basis_mo_value_device(qmckl_context context);

qmckl_compute_mo_basis_mo_value_device(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_value, double *restrict const mo_value);
