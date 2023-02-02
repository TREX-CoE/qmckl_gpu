#include "include/qmckl_ao.h"

//**********
// PROVIDE
//**********

/* shell_vgl */

qmckl_exit_code
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_provide_ao_basis_shell_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.shell_vgl_date) {

		/* Allocate array */
		if (ctx->ao_basis.shell_vgl == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->ao_basis.shell_num * 5 * ctx->point.num * sizeof(double);
			double *shell_vgl =
				(double *)qmckl_malloc_device(context, mem_info);

			if (shell_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_ao_basis_shell_vgl_device", NULL);
			}
			ctx->ao_basis.shell_vgl = shell_vgl;
		}

		qmckl_exit_code rc;
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
			return qmckl_failwith(context, QMCKL_FAILURE,
								  "compute_ao_basis_shell_vgl",
								  "Not yet implemented for basis type != 'G'");
		}
		if (rc != QMCKL_SUCCESS) {

			return rc;
		}

		ctx->ao_basis.shell_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

/* ao_vgl_gaussian */

qmckl_exit_code
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_ao_basis_ao_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

		qmckl_exit_code rc;

		/* Allocate array */
		if (ctx->ao_basis.ao_vgl == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double);
			double *ao_vgl = (double *)qmckl_malloc_device(context, mem_info);

			if (ao_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
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

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->ao_basis.ao_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

//**********
// GET
//**********

/* shell_vgl */

qmckl_exit_code
qmckl_get_ao_basis_shell_vgl_device(qmckl_context_device context,
									double *const shell_vgl,
									const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_vgl_device", NULL);
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_ao_basis_shell_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);
	// int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.shell_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_vgl",
							  "input array too small");
	}
	qmckl_memcpy_D2D(context, shell_vgl, ctx->ao_basis.shell_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS;
}

/* ao_vgl_gaussian */

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *const ao_vgl,
												 const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_ao_basis_ao_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_ao_vgl_device",
							  "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_vgl, ctx->ao_basis.ao_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS;
}
