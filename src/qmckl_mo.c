#include "include/qmckl_ao.h"

//**********
// PROVIDE
//**********

/* mo_vgl */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl_device(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check((qmckl_context) context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith((qmckl_context) context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_mo_basis_mo_vgl_device",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_mo_basis_mo_vgl_device",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->mo_basis.mo_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 5 * ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

    /* Allocate array */
    if (ctx->mo_basis.mo_vgl == NULL) {

      double* mo_vgl = (double*) qmckl_malloc_device(context, mem_info);

      if (mo_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_mo_basis_mo_vgl",
                               NULL);
      }
      ctx->mo_basis.mo_vgl = mo_vgl;
    }

    rc = qmckl_provide_ao_basis_ao_vgl_device(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_NOT_PROVIDED,
                             "qmckl_ao_basis",
                             NULL);
    }

    rc = qmckl_compute_mo_basis_mo_vgl_device(context,
                                       ctx->ao_basis.ao_num,
                                       ctx->mo_basis.mo_num,
                                       ctx->point.num,
                                       ctx->mo_basis.coefficient_t,
                                       ctx->ao_basis.ao_vgl,
                                       ctx->mo_basis.mo_vgl);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->mo_basis.mo_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

/* mo_value */

qmckl_exit_code qmckl_provide_mo_basis_mo_value_device(qmckl_context context) {

	qmckl_exit_code rc = QMCKL_SUCCESS;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_mo_basis_mo_value_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->mo_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_provide_mo_basis_mo_value_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->mo_basis.mo_value_date) {

		qmckl_exit_code rc;

		/* Allocate array */
		if (ctx->mo_basis.mo_value == NULL) {
			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

			double *mo_value = (double *)qmckl_malloc_device(context, mem_info);

			if (mo_value == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_mo_basis_mo_value", NULL);
			}
			ctx->mo_basis.mo_value = mo_value;
		}

		if (ctx->mo_basis.mo_vgl_date == ctx->point.date) {

			// mo_vgl has been computed at this step: Just copy the data.

			double *v = &(ctx->mo_basis.mo_value[0]);
			double *vgl = &(ctx->mo_basis.mo_vgl[0]);
#pragma omp target is_device_ptr(v, vgl)
			{
				for (int i = 0; i < ctx->point.num; ++i) {
					for (int k = 0; k < ctx->mo_basis.mo_num; ++k) {
						v[k] = vgl[k];
					}
					v += ctx->mo_basis.mo_num;
					vgl += ctx->mo_basis.mo_num * 5;
				}
			}

		} else {

			rc = qmckl_provide_ao_basis_ao_vgl_device(context);
			if (rc != QMCKL_SUCCESS) {
				return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
									  "qmckl_ao_basis_ao_vgl_device", NULL);
			}

			rc = qmckl_compute_mo_basis_mo_value_device(
				context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num,
				ctx->point.num, ctx->mo_basis.coefficient_t,
				ctx->ao_basis.ao_value, ctx->mo_basis.mo_value);
		}

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->mo_basis.mo_value_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

//**********
// GET
//**********

/* mo_vgl */

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_device(qmckl_context context,
												   double *const mo_vgl,
												   const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_mo_basis_mo_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = 5* ctx->point.num * ctx->mo_basis.mo_num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_vgl",
							  "input array too small");
	}
	qmckl_memcpy_D2D(context, mo_vgl, ctx->mo_basis.mo_vgl,
					 sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_inplace_device(
	qmckl_context context, double *const mo_vgl, const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_mo_basis_mo_vgl_device", NULL);
	}

	qmckl_exit_code rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = 5 * ctx->mo_basis.mo_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_vgl_device",
							  "input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	double *old_array = ctx->mo_basis.mo_vgl;

	ctx->mo_basis.mo_vgl = mo_vgl;

	rc = qmckl_provide_mo_basis_mo_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->mo_basis.mo_vgl = old_array;

	return QMCKL_SUCCESS;
}


/* mo_value */

qmckl_exit_code qmckl_get_mo_basis_mo_value_device(qmckl_context context,
												   double *const mo_value,
												   const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_mo_basis_mo_value_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->point.num * ctx->mo_basis.mo_num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_value",
							  "input array too small");
	}
	qmckl_memcpy_D2D(context, mo_value, ctx->mo_basis.mo_value,
					 sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_mo_basis_mo_value_inplace_device(
	qmckl_context context, double *const mo_value, const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_mo_basis_mo_value_device", NULL);
	}

	qmckl_exit_code rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->mo_basis.mo_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_value_device",
							  "input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	double *old_array = ctx->mo_basis.mo_value;

	ctx->mo_basis.mo_value = mo_value;

	rc = qmckl_provide_mo_basis_mo_value_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->mo_basis.mo_value = old_array;

	return QMCKL_SUCCESS;
}
