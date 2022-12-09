#pragma once

// This file contains prototypes for the AO functions

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

#include "qmckl_context_device.h"
#include "qmckl_memory_device.h"

qmckl_exit_code qmckl_ao_polynomial_transp_vgl_hpc_device(
	const qmckl_context_device context, const double *restrict X,
	const double *restrict R, const int32_t lmax, int64_t *restrict n,
	const int64_t ldl, double *restrict const VGL, const int64_t ldv);

qmckl_exit_code qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int32_t *restrict prim_num_per_nucleus,
	const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	const qmckl_matrix expo_per_nucleus, const qmckl_tensor coef_per_nucleus,
	double *restrict const ao_vgl);

qmckl_exit_code
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
									 double *const ao_vgl,
									 const int64_t size_max);
