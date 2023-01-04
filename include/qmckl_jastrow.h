#pragma once

// This file contains prototypes for the Jastrow functions
// TODO Get up to date with the latest Jastrow rework

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <qmckl.h>
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_jastrow_private_type.h"
#include "qmckl_jastrow_private_func.h"

qmckl_exit_code
qmckl_compute_tmp_c_offload(const qmckl_context context, const int64_t cord_num,
							const int64_t elec_num, const int64_t nucl_num,
							const int64_t walk_num,
							const double *een_rescaled_e,
							const double *een_rescaled_n, double *const tmp_c);

qmckl_exit_code qmckl_compute_dtmp_c_offload(
	const qmckl_context context, const int64_t cord_num, const int64_t elec_num,
	const int64_t nucl_num, const int64_t walk_num,
	const double *een_rescaled_e_deriv_e, const double *een_rescaled_n,
	double *const dtmp_c);

qmckl_exit_code qmckl_provide_tmp_c_offload(qmckl_context context);

qmckl_exit_code qmckl_provide_dtmp_c_offload(qmckl_context context);

qmckl_exit_code qmckl_provide_factor_een_deriv_e_offload(qmckl_context context);

qmckl_exit_code qmckl_get_jastrow_tmp_c_offload(qmckl_context context,
												double *const tmp_c);

qmckl_exit_code qmckl_get_jastrow_dtmp_c_offload(qmckl_context context,
												 double *const dtmp_c);

qmckl_exit_code
qmckl_get_jastrow_factor_een_deriv_e_offload(qmckl_context context,
											 double *const factor_een_deriv_e,
											 const int64_t size_max);
