#include "include/qmckl_ao.h"


qmckl_exit_code_device qmckl_init_ao_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_init_ao_basis_device", NULL);
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->ao_basis.uninitialized = (1 << 14) - 1;

	/* Default values */
	ctx->ao_basis.ao_cartesian = true;

	return QMCKL_SUCCESS_DEVICE;
}


//**********
// PROVIDE
//**********

/* shell_vgl */

qmckl_exit_code_device
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check_device((qmckl_context_device)context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
							  "qmckl_provide_ao_basis_shell_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.shell_vgl_date) {

		/* Allocate array */
		if (ctx->ao_basis.shell_vgl == NULL) {

			double *shell_vgl = (double *)qmckl_malloc_device(
				context,
				ctx->ao_basis.shell_num * 5 * ctx->point.num * sizeof(double));

			if (shell_vgl == NULL) {
				return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									  "qmckl_ao_basis_shell_vgl_device", NULL);
			}
			ctx->ao_basis.shell_vgl = shell_vgl;
		}

		qmckl_exit_code_device rc;
		if (ctx->ao_basis.type == 'G') {
			rc = qmckl_compute_ao_basis_shell_gaussian_vgl_device(
				context, ctx->ao_basis.prim_num, ctx->ao_basis.shell_num,
				ctx->point.num, ctx->nucleus.num,
				ctx->ao_basis.nucleus_shell_num, ctx->ao_basis.nucleus_index,
				ctx->ao_basis.nucleus_range, ctx->ao_basis.shell_prim_index,
				ctx->ao_basis.shell_prim_num, ctx->point.coord.data,
				ctx->nucleus.coord.data, ctx->ao_basis.exponent,
				ctx->ao_basis.coefficient_normalized, ctx->ao_basis.shell_vgl);
		} else {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
								  "compute_ao_basis_shell_vgl",
								  "Not yet implemented for basis type != 'G'");
		}
		if (rc != QMCKL_SUCCESS_DEVICE) {

			return rc;
		}

		ctx->ao_basis.shell_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_vgl */

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check_device((qmckl_context_device)context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith_device((qmckl_context_device)context, QMCKL_NOT_PROVIDED_DEVICE,
							  "qmckl_ao_basis_ao_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

		qmckl_exit_code_device rc;

		/* Allocate array */
		if (ctx->ao_basis.ao_vgl == NULL) {

			double *ao_vgl = (double *)qmckl_malloc_device(
				context,
				ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double));

			if (ao_vgl == NULL) {
				return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									  "qmckl_ao_basis_ao_vgl", NULL);
			}
			ctx->ao_basis.ao_vgl = ao_vgl;
		}

		/* Checking for shell_vgl */
		if (ctx->ao_basis.shell_vgl == NULL ||
			ctx->point.date > ctx->ao_basis.shell_vgl_date) {
			qmckl_provide_ao_basis_shell_vgl_device(context);
		}

		/* Compute ao_vgl_gaussian itself */
		if (ctx->ao_basis.type == 'G') {
			rc = qmckl_compute_ao_vgl_gaussian_device(
				context, ctx->ao_basis.ao_num, ctx->ao_basis.shell_num,
				ctx->point.num, ctx->nucleus.num, ctx->point.coord.data,
				ctx->nucleus.coord.data, ctx->ao_basis.nucleus_index,
				ctx->ao_basis.nucleus_shell_num, ctx->ao_basis.nucleus_range,
				ctx->ao_basis.nucleus_max_ang_mom, ctx->ao_basis.shell_ang_mom,
				ctx->ao_basis.ao_factor, ctx->ao_basis.shell_vgl,
				ctx->ao_basis.ao_vgl);
		} else {
			printf("Device pointers version of ao_vgl only "
				   "supports 'G' as its "
				   "ao_basis.type for now\n ");
		}

		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->ao_basis.ao_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GET
//**********

/* shell_vgl */

qmckl_exit_code_device
qmckl_get_ao_basis_shell_vgl_device(qmckl_context_device context,
									double *const shell_vgl,
									const int64_t size_max) {
	if (qmckl_context_check_device((qmckl_context_device)context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_get_ao_basis_shell_vgl_device", NULL);
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ao_basis_shell_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	// int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.shell_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_get_ao_basis_shell_vgl",
							  "input array too small");
	}
	qmckl_memcpy_D2D(context, shell_vgl, ctx->ao_basis.shell_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_vgl_gaussian */

qmckl_exit_code_device qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *const ao_vgl,
												 const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device((qmckl_context_device)context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_get_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_exit_code_device rc;
	rc = qmckl_provide_ao_basis_ao_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device((qmckl_context_device)context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_get_ao_basis_ao_vgl_device",
							  "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_vgl, ctx->ao_basis.ao_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_value */

qmckl_exit_code_device qmckl_get_ao_basis_ao_value_device(qmckl_context_device context,
												   double *const ao_value,
												   const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device((qmckl_context_device)context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_get_ao_basis_ao_value_device", NULL);
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ao_basis_ao_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device((qmckl_context_device)context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_get_ao_basis_ao_value_device",
							  "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_value, ctx->ao_basis.ao_value,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

/* Provided check  */

bool qmckl_ao_basis_provided(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->ao_basis.provided;
}


//**********
// VARIOUS GETTERS/SETTERS
//**********

qmckl_exit_code_device
qmckl_set_ao_basis_type_device(qmckl_context_device context, char basis_type) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_type_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (basis_type != 'G' && basis_type != 'S') {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_type_device", NULL);
	}

	int32_t mask = 1;

	ctx->ao_basis.type = basis_type;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (shell_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num <= 0");
	}

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (0L < prim_num && prim_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num > prim_num");
	}

	int32_t mask = 1 << 1;

	ctx->ao_basis.shell_num = shell_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t prim_num) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (prim_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "prim_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "shell_num is not set");
	}

	if (prim_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "prim_num < shell_num");
	}

	int32_t mask = 1 << 2;

	ctx->ao_basis.prim_num = prim_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_ao_num_device(qmckl_context_device context, int64_t ao_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_ao_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (ao_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "ao_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;
	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num is not set");
	}

	if (ao_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "ao_num < shell_num");
	}

	int32_t mask = 1 << 12;

	ctx->ao_basis.ao_num = ao_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 4;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 "nucl_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 "input array too small");
	}

	if (ctx->ao_basis.nucleus_index != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_index);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_nucleus_index_device", NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_index,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 3;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith_device(
			context, QMCKL_FAILURE_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device",
			"shell_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device",
			"input array too small");
	}

	if (ctx->ao_basis.nucleus_shell_num != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_shell_num);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_nucleus_shell_num_device",
				NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_ALLOCATION_FAILED_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_shell_num,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_shell_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 5;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num == 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_ang_mom != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_ang_mom);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
		}
	}

	int32_t *new_array =
		(int32_t *)qmckl_malloc_device(context, shell_num * sizeof(int32_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_ang_mom,
					 shell_num * sizeof(int32_t));

	ctx->ao_basis.shell_ang_mom = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_prim_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 6;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_prim_num != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_num);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_prim_num_device", NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_num,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 7;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(
			context, QMCKL_FAILURE_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device",
			"shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device",
			"input array too small");
	}

	if (ctx->ao_basis.shell_prim_index != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_index);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_prim_index_device",
				NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_ALLOCATION_FAILED_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_index,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 8;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, shell_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_factor,
					 shell_num * sizeof(double));

	ctx->ao_basis.shell_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 9;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 "input array too small");
	}

	if (ctx->ao_basis.exponent != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.exponent);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_exponent_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, exponent, prim_num * sizeof(double));

	ctx->ao_basis.exponent = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 10;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 "input array too small");
	}

	if (ctx->ao_basis.coefficient != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.coefficient);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_coefficient_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, coefficient,
					 prim_num * sizeof(double));

	ctx->ao_basis.coefficient = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	// Accepts an host array and copies it on device

	int32_t mask = 1 << 11;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.prim_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.prim_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_prim_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, prim_factor,
					 prim_num * sizeof(double));

	ctx->ao_basis.prim_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 13;

	int64_t ao_num = ctx->ao_basis.ao_num;

	if (ao_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 "ao_num is not set");
	}

	if (size_max < ao_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.ao_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.ao_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_ao_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, ao_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, ao_factor, ao_num * sizeof(double));

	ctx->ao_basis.ao_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *ao_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_ao_num", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 12;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_ao_num", NULL);
	}

	assert(ctx->ao_basis.ao_num > (int64_t)0);

	*ao_num = ctx->ao_basis.ao_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(const qmckl_context_device context,
									int64_t *const shell_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_factor", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_num", NULL);
	}

	assert(ctx->ao_basis.shell_num > (int64_t)0);
	*shell_num = ctx->ao_basis.shell_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(const qmckl_context_device context,
								   int64_t *prim_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_prim_num", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 2;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_prim_num", NULL);
	}

	assert(ctx->ao_basis.prim_num > (int64_t)0);

	*prim_num = ctx->ao_basis.prim_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_type_device(qmckl_context_device context,
							   char *const basis_type) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_type", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_type", NULL);
	}

	assert(ctx->ao_basis.type != (char)0);

	basis_type[0] = ctx->ao_basis.type;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 3;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(
			context, QMCKL_NOT_PROVIDED_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", NULL);
	}

	if (nucleus_shell_num == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_2_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device",
			"Array too small. Expected nucl_num");
	}

	assert(ctx->ao_basis.nucleus_shell_num != NULL);
	qmckl_memcpy_D2D(context, nucleus_shell_num,
					 ctx->ao_basis.nucleus_shell_num,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 4;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 NULL);
	}

	if (nucleus_index == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.nucleus_index != NULL);
	qmckl_memcpy_D2D(context, nucleus_index, ctx->ao_basis.nucleus_index,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 5;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	if (shell_ang_mom == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_ang_mom != NULL);
	qmckl_memcpy_D2D(context, shell_ang_mom, ctx->ao_basis.shell_ang_mom,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int32_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 8;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 NULL);
	}

	if (shell_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_factor != NULL);
	qmckl_memcpy_D2D(context, shell_factor, ctx->ao_basis.shell_factor,
					 (size_t)ctx->ao_basis.shell_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 6;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 NULL);
	}

	if (shell_prim_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_num != NULL);
	qmckl_memcpy_D2D(context, shell_prim_num, ctx->ao_basis.shell_prim_num,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 7;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(
			context, QMCKL_NOT_PROVIDED_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", NULL);
	}

	if (shell_prim_index == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_2_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device",
			"Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_index != NULL);
	qmckl_memcpy_D2D(context, shell_prim_index, ctx->ao_basis.shell_prim_index,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 9;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 NULL);
	}

	if (exponent == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.exponent != NULL);
	qmckl_memcpy_D2D(context, exponent, ctx->ao_basis.exponent,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 10;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 NULL);
	}

	if (coefficient == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 "Array too small. Expected prim_num");
	}
	assert(ctx->ao_basis.coefficient != NULL);
	qmckl_memcpy_D2D(context, coefficient, ctx->ao_basis.coefficient,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 11;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 NULL);
	}

	if (prim_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.prim_factor != NULL);
	qmckl_memcpy_D2D(context, prim_factor, ctx->ao_basis.prim_factor,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 13;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 NULL);
	}

	if (ao_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.ao_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 "Array too small. Expected ao_num");
	}

	assert(ctx->ao_basis.ao_factor != NULL);
	qmckl_memcpy_D2D(context, ao_factor, ctx->ao_basis.ao_factor,
					 ctx->ao_basis.ao_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}
