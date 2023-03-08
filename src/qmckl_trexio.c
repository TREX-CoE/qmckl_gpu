#include "../include/qmckl_trexio.h"

// This file provides wrappers to standard QMCkl functions accessible with the
// _device suffix. Only includes functions independent of OpenMP/OpenACC syntax.

//**********
// ELECTRON/POINT GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  int64_t up_num,
											  int64_t down_num) {
	return qmckl_set_electron_num((qmckl_context)context, up_num, down_num);
}

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context context,
												char transp, int64_t walk_num,
												double *coord,
												int64_t size_max) {

	size_t device_id = qmckl_get_device_id(context);
	int32_t mask = 0; // coord can be changed

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int64_t elec_num = ctx->electron.num;

	if (elec_num == 0L) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_electron_coord_device",
							  "elec_num is not set");
	}

	/* Swap pointers */
	qmckl_walker tmp = ctx->electron.walker_old;
	ctx->electron.walker_old = ctx->electron.walker;
	ctx->electron.walker = tmp;

	memcpy(&(ctx->point), &(ctx->electron.walker.point),
		   sizeof(qmckl_point_struct));

	qmckl_exit_code rc;
	rc = qmckl_set_point_device(context, transp, walk_num * elec_num, coord,
								size_max);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->electron.walker.num = walk_num;
	memcpy(&(ctx->electron.walker.point), &(ctx->point),
		   sizeof(qmckl_point_struct));

	return QMCKL_SUCCESS;
}

//**********
// NUCLEUS GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 int64_t num) {
	return qmckl_set_nucleus_num((qmckl_context)context, num);
}

qmckl_exit_code qmckl_get_nucleus_num_device(qmckl_context_device context,
											 int64_t *num) {
	return qmckl_get_nucleus_num((qmckl_context)context, num);
}

qmckl_exit_code qmckl_set_nucleus_charge_device(qmckl_context_device context,
												double *charge,
												int64_t size_max) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_charge_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (charge == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_nucleus_charge_device",
							  "charge is a null pointer");
	}

	int64_t num;
	qmckl_exit_code rc;

	rc = qmckl_get_nucleus_num(context, &num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	if (num > size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_nucleus_charge_device",
							  "Array too small");
	}

	ctx->nucleus.charge = qmckl_vector_alloc_device(context, num);
	rc = qmckl_vector_of_double_device(context, charge, num,
									   &(ctx->nucleus.charge));

	if (rc != QMCKL_SUCCESS) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_nucleus_charge_device",
							  "Error in vector->double* conversion");
	}

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

//**********
// AO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   char basis_type) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_type_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (basis_type != 'G' && basis_type != 'S') {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_type_device", NULL);
	}

	int32_t mask = 1;

	ctx->ao_basis.type = basis_type;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_num_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (shell_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num <= 0");
	}

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (0L < prim_num && prim_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num > prim_num");
	}

	int32_t mask = 1 << 1;

	ctx->ao_basis.shell_num = shell_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t prim_num) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_num_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (prim_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_prim_num_device",
							  "prim_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_prim_num_device",
							  "shell_num is not set");
	}

	if (prim_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_prim_num_device",
							  "prim_num < shell_num");
	}

	int32_t mask = 1 << 2;

	ctx->ao_basis.prim_num = prim_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t ao_num) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_ao_num_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (ao_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "ao_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;
	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num is not set");
	}

	if (ao_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "ao_num < shell_num");
	}

	int32_t mask = 1 << 12;

	ctx->ao_basis.ao_num = ao_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_nucleus_index_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 4;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_index_device",
							  "nucl_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_index_device",
							  "input array too small");
	}

	if (ctx->ao_basis.nucleus_index != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_index);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_nucleus_index_device",
								  NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_nucleus_index_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_index,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 3;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  "shell_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  "input array too small");
	}

	if (ctx->ao_basis.nucleus_shell_num != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_shell_num);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith(context, rc,
								  "qmckl_set_ao_basis_nucleus_shell_num_device",
								  NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_shell_num,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_shell_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 5;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num == 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_ang_mom_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_ang_mom_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_ang_mom != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_ang_mom);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_ang_mom_device",
								  NULL);
		}
	}

	int32_t *new_array =
		(int32_t *)qmckl_malloc_device(context, shell_num * sizeof(int32_t));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_ang_mom,
					 shell_num * sizeof(int32_t));

	ctx->ao_basis.shell_ang_mom = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(
			(qmckl_context)context, QMCKL_INVALID_CONTEXT,
			"qmckl_set_ao_basis_nucleus_shell_prim_num_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 6;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_prim_num_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_prim_num_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_prim_num != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_num);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_prim_num_device",
								  NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_prim_num_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_num,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 7;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_prim_index != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_index);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_prim_index_device",
								  NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_index,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_factor_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 8;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_factor_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_factor_device",
								  NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, shell_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_factor_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_factor,
					 shell_num * sizeof(double));

	ctx->ao_basis.shell_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_exponent_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 9;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_exponent_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_exponent_device",
							  "input array too small");
	}

	if (ctx->ao_basis.exponent != NULL) {
		qmckl_exit_code rc = qmckl_free_device(context, ctx->ao_basis.exponent);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_exponent_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_exponent_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, exponent, prim_num * sizeof(double));

	ctx->ao_basis.exponent = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_coefficient_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 10;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_coefficient_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_coefficient_device",
							  "input array too small");
	}

	if (ctx->ao_basis.coefficient != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.coefficient);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_coefficient_device",
								  NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_coefficient_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, coefficient,
					 prim_num * sizeof(double));

	ctx->ao_basis.coefficient = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_prim_factor_device(qmckl_context context,
													  double *prim_factor,
													  int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	// Accepts an host array and copies it on device

	int32_t mask = 1 << 11;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_prim_factor_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_prim_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.prim_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.prim_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_prim_factor_device",
								  NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, prim_factor,
					 prim_num * sizeof(double));

	ctx->ao_basis.prim_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}
	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 13;

	int64_t ao_num = ctx->ao_basis.ao_num;

	if (ao_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_ao_factor_device",
							  "ao_num is not set");
	}

	if (size_max < ao_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_ao_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.ao_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.ao_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_ao_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, ao_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_ao_factor_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, ao_factor, ao_num * sizeof(double));

	ctx->ao_basis.ao_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_ao_num_device(qmckl_context_device context,
												 int64_t *ao_num) {
	return qmckl_get_ao_basis_ao_num((qmckl_context)context, ao_num);
}

qmckl_exit_code
qmckl_get_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t *ao_num) {
	return qmckl_get_ao_basis_shell_num((qmckl_context)context, ao_num);
}

qmckl_exit_code qmckl_get_ao_basis_prim_num_device(qmckl_context_device context,
												   int64_t *prim_num) {
	return qmckl_get_ao_basis_prim_num((qmckl_context)context, prim_num);
}

qmckl_exit_code qmckl_get_ao_basis_type_device(qmckl_context_device context,
											   char *type) {
	return qmckl_get_ao_basis_type((qmckl_context)context, type);
}

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_nucleus_shell_num_device",
							  NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 3;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_nucleus_shell_num_device",
							  NULL);
	}

	if (nucleus_shell_num == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_nucleus_shell_num_device",
							  "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_nucleus_shell_num_device",
							  "Array too small. Expected nucl_num");
	}

	assert(ctx->ao_basis.nucleus_shell_num != NULL);
	qmckl_memcpy_D2D(context, nucleus_shell_num,
					 ctx->ao_basis.nucleus_shell_num,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_nucleus_index", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 4;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_nucleus_index_device", NULL);
	}

	if (nucleus_index == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_nucleus_index_device",
							  "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_nucleus_index_device",
							  "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.nucleus_index != NULL);
	qmckl_memcpy_D2D(context, nucleus_index, ctx->ao_basis.nucleus_index,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_ang_mom_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 5;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_shell_ang_mom_device", NULL);
	}

	if (shell_ang_mom == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_shell_ang_mom_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_ang_mom_device",
							  "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_ang_mom != NULL);
	qmckl_memcpy_D2D(context, shell_ang_mom, ctx->ao_basis.shell_ang_mom,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int32_t));
	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_factor_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 8;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_shell_factor_device", NULL);
	}

	if (shell_factor == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_shell_factor_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_factor_device",
							  "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_factor != NULL);
	qmckl_memcpy_D2D(context, shell_factor, ctx->ao_basis.shell_factor,
					 (size_t)ctx->ao_basis.shell_num * sizeof(double));
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, double *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_prim_num_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 6;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_shell_prim_num_device", NULL);
	}

	if (shell_prim_num == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_shell_prim_num_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_prim_num_device",
							  "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_num != NULL);
	qmckl_memcpy_D2D(context, shell_prim_num, ctx->ao_basis.shell_prim_num,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_prim_index_device",
							  NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 7;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_shell_prim_index_device",
							  NULL);
	}

	if (shell_prim_index == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_shell_prim_index_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_prim_index_device",
							  "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_index != NULL);
	qmckl_memcpy_D2D(context, shell_prim_index, ctx->ao_basis.shell_prim_index,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
												   double *exponent,
												   int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_exponent_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 9;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_exponent_device", NULL);
	}

	if (exponent == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_exponent_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_exponent_device",
							  "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.exponent != NULL);
	qmckl_memcpy_D2D(context, exponent, ctx->ao_basis.exponent,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_coefficient_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 10;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_coefficient_device", NULL);
	}

	if (coefficient == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_coefficient_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_coefficient_device",
							  "Array too small. Expected prim_num");
	}
	assert(ctx->ao_basis.coefficient != NULL);
	qmckl_memcpy_D2D(context, coefficient, ctx->ao_basis.coefficient,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_prim_factor_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 11;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_prim_factor_device", NULL);
	}

	if (prim_factor == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_prim_factor_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_prim_factor_device",
							  "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.prim_factor != NULL);
	qmckl_memcpy_D2D(context, prim_factor, ctx->ao_basis.prim_factor,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_factor_device", NULL);
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 13;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_ao_factor_device", NULL);
	}

	if (ao_factor == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_ao_basis_ao_factor_device",
							  "NULL pointer");
	}

	if (size_max < ctx->ao_basis.ao_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_ao_factor_device",
							  "Array too small. Expected ao_num");
	}

	assert(ctx->ao_basis.ao_factor != NULL);
	qmckl_memcpy_D2D(context, ao_factor, ctx->ao_basis.ao_factor,
					 ctx->ao_basis.ao_num * sizeof(double));
	return QMCKL_SUCCESS;
}

//**********
// MO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_mo_basis_mo_num_device(qmckl_context_device context,
												 int64_t mo_num) {

	int32_t mask = 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALREADY_SET,
							  "qmckl_set_mo_basis_mo_num_device", NULL);
	}

	if (mo_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_mo_basis_mo_num_device",
							  "mo_num <= 0");
	}

	ctx->mo_basis.mo_num = mo_num;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_num_device(const qmckl_context_device context,
								 int64_t *mo_num) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_mo_basis_mo_num_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->mo_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_mo_basis_mo_num", NULL);
	}

	assert(ctx->mo_basis.mo_num > (int64_t)0);
	*mo_num = ctx->mo_basis.mo_num;
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_mo_basis_coefficient_device(qmckl_context context,
													  double *coefficient) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALREADY_SET,
							  "qmckl_set_mo_basis_coefficient_device", NULL);
	}

	if (ctx->mo_basis.coefficient != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->mo_basis.coefficient);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_mo_basis_coefficient_device",
								  NULL);
		}
	}

	double *new_array = (double *)qmckl_malloc_device(
		context, ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double));
	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_mo_basis_coefficient_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, coefficient,
					 ctx->ao_basis.ao_num * ctx->mo_basis.mo_num *
						 sizeof(double));

	ctx->mo_basis.coefficient = new_array;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

//**********
// TREXIO CONTEXT FILL
//**********

qmckl_exit_code
qmckl_trexio_read_electron_X_device(qmckl_context_device context,
									trexio_t *file) {

	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	int rcio = 0;

	int64_t up_num = 0L;
	int64_t dn_num = 0L;

	rcio = trexio_read_electron_up_num_64(file, &up_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_electron_up_num",
							  trexio_string_of_error(rcio));
	}

	assert(up_num >= 0L);

	rcio = trexio_read_electron_dn_num_64(file, &dn_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_electron_dn_num",
							  trexio_string_of_error(rcio));
	}

	assert(dn_num >= 0L);

	qmckl_exit_code rc;
	rc = qmckl_set_electron_num_device(context, up_num, dn_num);
	return rc;
}

qmckl_exit_code qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
												   trexio_t *file) {
	assert(context != (qmckl_context)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;

	int64_t nucleus_num = 0L;

	rcio = trexio_read_nucleus_num_64(file, &nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_nucleus_num",
							  trexio_string_of_error(rcio));
	}

	assert(nucleus_num > 0);
	rc = qmckl_set_nucleus_num_device(context, nucleus_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = nucleus_num * sizeof(double);

		double *nucl_charge_h = (double *)qmckl_malloc_host(context, mem_info);
		double *nucl_charge_d = (double *)qmckl_malloc_device(
			context, nucleus_num * sizeof(double));

		if (nucl_charge_h == NULL || nucl_charge_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_nucleus_X_device", NULL);
		}

		assert(nucl_charge_h != NULL && nucl_charge_d != NULL);

		rcio = trexio_read_safe_nucleus_charge_64(file, nucl_charge_h,
												  nucleus_num);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_nucleus_charge",
								  trexio_string_of_error(rcio));
		}

		qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge_h, mem_info.size);
		rc = qmckl_set_nucleus_charge_device(context, nucl_charge_d,
											 nucleus_num);

		qmckl_free_host(context, nucl_charge_h);
		qmckl_free_device(context, nucl_charge_d);
		nucl_charge_h = NULL;
		nucl_charge_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = nucleus_num * 3 * sizeof(double);

	double *nucl_coord_h = (double *)qmckl_malloc_host(context, mem_info);
	double *nucl_coord_d = (double *)qmckl_malloc_device(
		context, nucleus_num * 3 * sizeof(double));

	if (nucl_coord_h == NULL || nucl_coord_d == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_trexio_read_nucleus_X_device", NULL);
	}

	assert(nucl_coord_h != NULL && nucl_coord_d != NULL);

	rcio =
		trexio_read_safe_nucleus_coord_64(file, nucl_coord_h, 3 * nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_nucleus_charge",
							  trexio_string_of_error(rcio));
	}

	qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord_h, mem_info.size);
	rc = qmckl_set_nucleus_coord_device(context, 'N', nucl_coord_d,
										3 * nucleus_num);

	qmckl_free_host(context, nucl_coord_h);
	qmckl_free_device(context, nucl_coord_d);
	nucl_coord_h = NULL;
	nucl_coord_d = NULL;

	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_ao_X_device(qmckl_context context,
											  trexio_t *file) {
	assert(context != (qmckl_context)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;
	int64_t nucleus_num = 0L;

	rc = qmckl_get_nucleus_num_device(context, &nucleus_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

#define MAX_STR_LEN 1024
	char basis_type[MAX_STR_LEN];

	rcio = trexio_read_basis_type(file, basis_type, MAX_STR_LEN);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_type",
							  trexio_string_of_error(rcio));
	}

	if (basis_type[0] == 'G') {
		rc = qmckl_set_ao_basis_type_device(context, basis_type[0]);
	} else {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_type", "Invalid basis type");
	}

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t shell_num = 0L;

	rcio = trexio_read_basis_shell_num_64(file, &shell_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_shell_num",
							  trexio_string_of_error(rcio));
	}

	assert(shell_num > 0);
	rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t prim_num = 0L;

	rcio = trexio_read_basis_prim_num_64(file, &prim_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_prim_num",
							  trexio_string_of_error(rcio));
	}

	assert(prim_num > 0);
	rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t ao_num = 0LL;

	rcio = trexio_read_ao_num_64(file, &ao_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_ao_num",
							  trexio_string_of_error(rcio));
	}

	assert(ao_num > 0);
	rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;
		int64_t *nucleus_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_index_d = (int64_t *)qmckl_malloc_device(
			context, nucleus_num * sizeof(int64_t));

		if (nucleus_index_h == NULL || nucleus_index_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "index_X_device",
								  NULL);
		}

		assert(nucleus_index_h != NULL && nucleus_index_d != NULL);

		/* Allocate temporary array */
		mem_info.size = shell_num * sizeof(int64_t);
		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "index_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read in the temporary array */
		rcio =
			trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_nucleus_index",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
													 nucleus_num);
		if (rc != QMCKL_SUCCESS) {
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return rc;
		}

		for (int i = shell_num - 1; i >= 0; --i) {
			int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_index_h);
				qmckl_free_device(context, nucleus_index_d);
				nucleus_index_h = NULL;
				nucleus_index_d = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_nucleus_index",
									  "Irrelevant data in TREXIO file");
			}
			nucleus_index_h[k] = i;
		}

		qmckl_memcpy_H2D(context, nucleus_index_d, nucleus_index_h,
						 size_backup);

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
													 shell_num);

		qmckl_free_host(context, nucleus_index_h);
		qmckl_free_device(context, nucleus_index_d);
		nucleus_index_h = NULL;
		nucleus_index_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *nucleus_shell_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_shell_num_d = (int64_t *)qmckl_malloc_device(
			context, nucleus_num * sizeof(int64_t));

		if (nucleus_shell_num_h == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "shell_num_X_device",
								  NULL);
		}

		assert(nucleus_shell_num_h != NULL);

		/* Allocate temporary array */
		mem_info.size = shell_num * sizeof(int64_t);
		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, nucleus_shell_num_h);
			qmckl_free_device(context, nucleus_shell_num_d);
			nucleus_shell_num_h = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "shell_num_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read in the temporary array */
		rcio =
			trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			qmckl_free_host(context, nucleus_shell_num_h);
			qmckl_free_device(context, nucleus_shell_num_d);
			nucleus_shell_num_h = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_nucleus_shell_num",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < nucleus_num; ++i) {
			nucleus_shell_num_h[i] = 0;
		}

		for (int i = 0; i < shell_num; ++i) {
			int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_shell_num_h);
				qmckl_free_device(context, nucleus_shell_num_d);
				nucleus_shell_num_h = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_nucleus_shell_num",
									  "Irrelevant data in TREXIO file");
			}
			nucleus_shell_num_h[k] += 1;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, nucleus_shell_num_d, nucleus_shell_num_h,
						 size_backup);
		rc = qmckl_set_ao_basis_nucleus_shell_num_device(
			context, nucleus_shell_num_d, shell_num);

		qmckl_free_host(context, nucleus_shell_num_h);
		qmckl_free_device(context, nucleus_shell_num_d);
		nucleus_shell_num_h = NULL;
		nucleus_shell_num_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int32_t);

		int32_t *shell_ang_mom_h =
			(int32_t *)qmckl_malloc_host(context, mem_info);
		int32_t *shell_ang_mom_d = (int32_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int32_t));

		if (shell_ang_mom_h == NULL || shell_ang_mom_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "ang_mom_X_device",
								  NULL);
		}

		assert(shell_ang_mom_h != NULL && shell_ang_mom_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_ang_mom_32(file, shell_ang_mom_h,
													   shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_ang_mom_h);
			qmckl_free_device(context, shell_ang_mom_d);
			shell_ang_mom_h = NULL;
			shell_ang_mom_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_ang_mom",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, shell_ang_mom_d, shell_ang_mom_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom_d,
													 shell_num);

		qmckl_free_host(context, shell_ang_mom_h);
		qmckl_free_device(context, shell_ang_mom_d);
		shell_ang_mom_h = NULL;
		shell_ang_mom_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_num_d = (int64_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int64_t));

		if (shell_prim_num_h == NULL || shell_prim_num_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_num_X_device",
								  NULL);
		}

		assert(shell_prim_num_h != NULL && shell_prim_num_d != NULL);

		/* Allocate temporary array */
		mem_info.size = prim_num * sizeof(int64_t);

		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, shell_prim_num_h);
			qmckl_free_device(context, shell_prim_num_d);
			shell_prim_num_h = NULL;
			shell_prim_num_d = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_num_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_prim_num_h);
			qmckl_free_device(context, shell_prim_num_d);
			shell_prim_num_h = NULL;
			shell_prim_num_d = NULL;
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_prim_num",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < shell_num; ++i) {
			shell_prim_num_h[i] = 0;
		}

		for (int i = 0; i < prim_num; ++i) {
			int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				qmckl_free_host(context, shell_prim_num_h);
				qmckl_free_device(context, shell_prim_num_d);
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_shell_prim_num",
									  "Irrelevant data in TREXIO file");
			}
			shell_prim_num_h[k] += 1;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, shell_prim_num_d, shell_prim_num_h,
						 size_backup);
		rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num_d,
													  shell_num);

		qmckl_free_host(context, shell_prim_num_h);
		qmckl_free_device(context, shell_prim_num_d);
		shell_prim_num_h = NULL;
		shell_prim_num_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_index_d = (int64_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int64_t));

		if (shell_prim_index_h == NULL || shell_prim_index_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_index_X_device",
								  NULL);
		}

		assert(shell_prim_index_h != NULL && shell_prim_index_d != NULL);

		/* Allocate temporary array */
		mem_info.size = prim_num * sizeof(int64_t);

		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_index_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_prim_index_h);
			qmckl_free_device(context, shell_prim_index_d);
			shell_prim_index_h = NULL;
			shell_prim_index_d = NULL;
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_prim_index",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = prim_num - 1; i >= 0; --i) {
			int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, shell_prim_index_h);
				qmckl_free_device(context, shell_prim_index_d);
				shell_prim_index_h = NULL;
				shell_prim_index_d = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_shell_prim_index",
									  "Irrelevant data in TREXIO file");
			}
			shell_prim_index_h[k] = i;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, shell_prim_index_d, shell_prim_index_h,
						 size_backup);
		rc = qmckl_set_ao_basis_shell_prim_index_device(
			context, shell_prim_index_d, shell_num);

		qmckl_free_host(context, shell_prim_index_h);
		qmckl_free_device(context, shell_prim_index_d);
		shell_prim_index_h = NULL;
		shell_prim_index_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(double);

		double *shell_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *shell_factor_d =
			(double *)qmckl_malloc_device(context, shell_num * sizeof(double));

		if (shell_factor_h == NULL || shell_factor_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_shell_factor_X_device", NULL);
		}

		assert(shell_factor_h != NULL && shell_factor_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_factor_64(file, shell_factor_h,
													  shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_factor_h);
			qmckl_free_device(context, shell_factor_d);
			shell_factor_h = NULL;
			shell_factor_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_factor",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, shell_factor_d, shell_factor_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor_d,
													shell_num);

		qmckl_free_host(context, shell_factor_h);
		qmckl_free_device(context, shell_factor_d);
		shell_factor_h = NULL;
		shell_factor_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *exponent_h = (double *)qmckl_malloc_host(context, mem_info);
		double *exponent_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (exponent_h == NULL || exponent_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_exponent_X", NULL);
		}

		assert(exponent_h != NULL && exponent_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_exponent_64(file, exponent_h, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, exponent_h);
			qmckl_free_device(context, exponent_d);
			exponent_h = NULL;
			exponent_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_exponent",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, exponent_d, exponent_h, mem_info.size);
		rc = qmckl_set_ao_basis_exponent_device(context, exponent_d, prim_num);

		qmckl_free_host(context, exponent_h);
		qmckl_free_device(context, exponent_d);
		exponent_h = NULL;
		exponent_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *coefficient_h = (double *)qmckl_malloc_host(context, mem_info);
		double *coefficient_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (coefficient_h == NULL || coefficient_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_coefficient_X_device", NULL);
		}

		assert(coefficient_h != NULL && coefficient_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_coefficient_64(file, coefficient_h,
													 prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, coefficient_h);
			qmckl_free_device(context, coefficient_d);
			coefficient_h = NULL;
			coefficient_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_coefficient",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, coefficient_d, coefficient_h, mem_info.size);
		rc = qmckl_set_ao_basis_coefficient_device(context, coefficient_d,
												   prim_num);

		qmckl_free_host(context, coefficient_h);
		qmckl_free_device(context, coefficient_d);
		coefficient_h = NULL;
		coefficient_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *prim_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *prim_factor_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (prim_factor_h == NULL || prim_factor_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_prim_factor_X_device", NULL);
		}

		assert(prim_factor_h != NULL && prim_factor_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_prim_factor_64(file, prim_factor_h,
													 prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, prim_factor_h);
			qmckl_free_device(context, prim_factor_d);
			prim_factor_h = NULL;
			prim_factor_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_prim_factor",
								  trexio_string_of_error(rcio));
		}

		/* Read data */
		qmckl_memcpy_H2D(context, prim_factor_d, prim_factor_h, mem_info.size);
		rc = qmckl_set_ao_basis_prim_factor_device(context, prim_factor_d,
												   prim_num);

		qmckl_free_host(context, prim_factor_h);
		qmckl_free_device(context, prim_factor_d);
		prim_factor_h = NULL;
		prim_factor_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = ao_num * sizeof(double);

		double *ao_normalization_h =
			(double *)qmckl_malloc_host(context, mem_info);
		double *ao_normalization_d =
			(double *)qmckl_malloc_device(context, ao_num * sizeof(double));

		if (ao_normalization_h == NULL || ao_normalization_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_ao_normalization_X_device", NULL);
		}

		assert(ao_normalization_h != NULL && ao_normalization_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_ao_normalization_64(file, ao_normalization_h,
													ao_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, ao_normalization_h);
			qmckl_free_device(context, ao_normalization_d);
			ao_normalization_h = NULL;
			ao_normalization_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_ao_normalization",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, ao_normalization_d, ao_normalization_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_ao_factor_device(context, ao_normalization_d,
												 ao_num);

		qmckl_free_host(context, ao_normalization_h);
		qmckl_free_device(context, ao_normalization_d);
		ao_normalization_h = NULL;
		ao_normalization_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_mo_X_device(qmckl_context_device context,
											  trexio_t *file) {
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;
	int64_t ao_num = 0L;

	rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t mo_num = 0L;

	rcio = trexio_read_mo_num_64(file, &mo_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_mo_num",
							  trexio_string_of_error(rcio));
	}

	assert(mo_num > 0);
	rc = qmckl_set_mo_basis_mo_num_device(context, mo_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ao_num * mo_num * sizeof(double);

		double *mo_coef_h = (double *)qmckl_malloc_host(context, mem_info);
		double *mo_coef_d = (double *)qmckl_malloc_device(
			context, ao_num * mo_num * sizeof(double));

		if (mo_coef_h == NULL || mo_coef_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_mo_X_device", NULL);
		}

		assert(mo_coef_h != NULL && mo_coef_d != NULL);

		rcio = trexio_read_mo_coefficient_64(file, mo_coef_h);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_mo_coefficient",
								  trexio_string_of_error(rcio));
		}

		qmckl_memcpy_H2D(context, mo_coef_d, mo_coef_h, mem_info.size);
		rc = qmckl_set_mo_basis_coefficient_device(context, mo_coef_d);

		qmckl_free_host(context, mo_coef_h);
		qmckl_free_device(context, mo_coef_d);
		mo_coef_h = NULL;
		mo_coef_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_device(qmckl_context_device context,
										 char *file_name, int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_exit_code rc;
	char file_name_new[size_max + 1];
	strncpy(file_name_new, file_name, size_max + 1);
	file_name_new[size_max] = '\0';

	trexio_t *file = qmckl_trexio_open_X(file_name_new, &rc);
	if (file == NULL) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_trexio_read_device",
							  trexio_string_of_error(rc));
	}

	assert(file != NULL);

	rc = qmckl_trexio_read_electron_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device",
							  "Error reading electron");
	}

	rc = qmckl_trexio_read_nucleus_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device",
							  "Error reading nucleus");
	}

	rc = qmckl_trexio_read_ao_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device", "Error reading AOs");
	}

	rc = qmckl_trexio_read_mo_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_omp_read", "Error reading MOs");
	}

	trexio_close(file);
	file = NULL;
	return rc;
}
