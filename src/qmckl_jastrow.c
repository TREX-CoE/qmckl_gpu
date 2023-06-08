#include "../include/qmckl_jastrow.h"

/* Initialized check */

bool qmckl_jastrow_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->jastrow.provided;
}

//**********
// FINALIZE
//**********

qmckl_exit_code_device
qmckl_finalize_jastrow_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* ----------------------------------- */
	/* Check for the necessary information */
	/* ----------------------------------- */

	if (!(ctx->electron.provided)) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_electron", NULL);
	}

	if (!(ctx->nucleus.provided)) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_nucleus", NULL);
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_asymp_jasa_device(context);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_provide_jastrow_asymp_jasb_device(context);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_context_touch_device(context);
	return rc;
}

//**********
// FINALIZE
//**********

qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled) {

	int k;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	for (int k = 0; k < walk_num; k++) {
		info = qmckl_distance_rescaled_device(
			context, 'T', 'T', elec_num, elec_num, coord + (k * elec_num),
			elec_num * walk_num, coord + (k * elec_num), elec_num * walk_num,
			ee_distance_rescaled + (k * elec_num * elec_num), elec_num,
			rescale_factor_ee);
		if (info != QMCKL_SUCCESS_DEVICE) {
			break;
		}
	}
	return info;
}

qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_deriv_e_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled_deriv_e) {

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	for (int k = 0; k < walk_num; k++) {
		info = qmckl_distance_rescaled_deriv_e_device(
			context, 'T', 'T', elec_num, elec_num, coord + (k * elec_num),
			elec_num * walk_num, coord + (k * elec_num), elec_num * walk_num,
			ee_distance_rescaled_deriv_e + (k * 4 * elec_num * elec_num),
			elec_num, rescale_factor_ee);

		if (info != QMCKL_SUCCESS_DEVICE)
			break;
	}

	return info;
}

//**********
// SETTERS
//**********

qmckl_exit_code_device
qmckl_set_jastrow_rescale_factor_ee_device(qmckl_context_device context,
										   const double rescale_factor_ee) {
	int32_t mask = 1 << 8;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	if (rescale_factor_ee <= 0.0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_rescale_factor_ee",
									 "rescale_factor_ee <= 0.0");
	}

	ctx->jastrow.rescale_factor_ee = rescale_factor_ee;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_aord_num_device(qmckl_context_device context,
								  const int64_t aord_num) {
	int32_t mask = 1 << 0;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	if (aord_num < 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_aord_num",
									 "aord_num < 0");
	}

	ctx->jastrow.aord_num = aord_num;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_bord_num_device(qmckl_context_device context,
								  const int64_t bord_num) {
	int32_t mask = 1 << 1;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	if (bord_num < 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_bord_num",
									 "bord_num < 0");
	}

	ctx->jastrow.bord_num = bord_num;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_cord_num_device(qmckl_context_device context,
								  const int64_t cord_num) {
	int32_t mask = 1 << 2;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	if (cord_num < 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_cord_num",
									 "cord_num < 0");
	}

	int64_t dim_c_vector = -1;
	qmckl_exit_code_device rc =
		qmckl_compute_dim_c_vector_device(context, cord_num, &dim_c_vector);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	ctx->jastrow.cord_num = cord_num;
	ctx->jastrow.dim_c_vector = dim_c_vector;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_type_nucl_num_device(qmckl_context_device context,
									   const int64_t type_nucl_num) {
	int32_t mask = 1 << 3;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	if (type_nucl_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_type_nucl_num",
									 "type_nucl_num < 0");
	}

	ctx->jastrow.type_nucl_num = type_nucl_num;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_type_nucl_vector_device(qmckl_context_device context,
										  const int64_t *type_nucl_vector,
										  const int64_t nucl_num) {
	int32_t mask = 1 << 4;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	int64_t type_nucl_num = ctx->jastrow.type_nucl_num;

	if (type_nucl_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_type_nucl_vector",
									 "type_nucl_num not initialized");
	}

	if (type_nucl_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_type_nucl_vector",
									 "type_nucl_vector = NULL");
	}

	if (ctx->jastrow.type_nucl_vector != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->jastrow.type_nucl_vector);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_type_nucl_vector",
				"Unable to free ctx->jastrow.type_nucl_vector");
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_jastrow_type_nucl_vector",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, type_nucl_vector,
					 nucl_num * sizeof(int64_t));

	ctx->jastrow.type_nucl_vector = new_array;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_a_vector_device(qmckl_context_device context,
								  const double *a_vector,
								  const int64_t size_max) {
	int32_t mask = 1 << 5;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	int64_t aord_num = ctx->jastrow.aord_num;
	if (aord_num < 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_a_vector",
									 "aord_num not initialized");
	}

	int64_t type_nucl_num = ctx->jastrow.type_nucl_num;

	if (type_nucl_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_a_vector",
									 "type_nucl_num not initialized");
	}

	if (a_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_a_vector",
									 "a_vector = NULL");
	}

	if (ctx->jastrow.a_vector != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->jastrow.a_vector);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_jastrow_a_vector",
				"Unable to free ctx->jastrow.a_vector");
		}
	}

	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;
	mem_info.size = (aord_num + 1) * type_nucl_num * sizeof(double);

	if (size_max < (aord_num + 1) * type_nucl_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE, "qmckl_set_jastrow_a_vector",
			"Array too small. Expected (aord_num+1)*type_nucl_num");
	}

	double *new_array = (double *)qmckl_malloc_device(context, mem_info.size);

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_jastrow_coefficient", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, a_vector,
					 (aord_num + 1) * type_nucl_num * sizeof(double));

	ctx->jastrow.a_vector = new_array;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_b_vector_device(qmckl_context_device context,
								  const double *b_vector,
								  const int64_t size_max) {
	int32_t mask = 1 << 6;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	int64_t bord_num = ctx->jastrow.bord_num;
	if (bord_num < 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_b_vector",
									 "bord_num not initialized");
	}

	if (b_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_b_vector",
									 "b_vector = NULL");
	}

	if (ctx->jastrow.b_vector != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->jastrow.b_vector);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_jastrow_b_vector",
				"Unable to free ctx->jastrow.b_vector");
		}
	}

	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;
	mem_info.size = (bord_num + 1) * sizeof(double);

	if (size_max < (bord_num + 1)) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_jastrow_b_vector",
									 "Array too small. Expected (bord_num+1)");
	}

	double *new_array = (double *)qmckl_malloc_device(context, mem_info.size);

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_jastrow_coefficient", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, b_vector,
					 (bord_num + 1) * sizeof(double));

	ctx->jastrow.b_vector = new_array;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_jastrow_c_vector_device(qmckl_context_device context,
								  const double *c_vector,
								  const int64_t size_max) {
	int32_t mask = 1 << 7;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_jastrow_*", NULL);
	}

	int64_t type_nucl_num = ctx->jastrow.type_nucl_num;
	if (type_nucl_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_c_vector",
									 "type_nucl_num not initialized");
	}

	int64_t dim_c_vector = ctx->jastrow.dim_c_vector;
	if (dim_c_vector < 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_set_jastrow_c_vector",
									 "cord_num not initialized");
	}

	if (c_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_jastrow_c_vector",
									 "c_vector = NULL");
	}

	if (ctx->jastrow.c_vector != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->jastrow.c_vector);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_jastrow_c_vector",
				"Unable to free ctx->jastrow.c_vector");
		}
	}

	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;
	mem_info.size = dim_c_vector * type_nucl_num * sizeof(double);

	if ((size_t)size_max < dim_c_vector * type_nucl_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE, "qmckl_set_jastrow_c_vector",
			"Array too small. Expected dim_c_vector * type_nucl_num");
	}

	double *new_array = (double *)qmckl_malloc_device(context, mem_info.size);

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_jastrow_coefficient", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, c_vector,
					 dim_c_vector * type_nucl_num * sizeof(double));

	ctx->jastrow.c_vector = new_array;

	ctx->jastrow.uninitialized &= ~mask;
	ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
	if (ctx->jastrow.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_jastrow_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GETTERS (basic)
//**********

qmckl_exit_code_device
qmckl_get_jastrow_aord_num_device(qmckl_context_device context,
								  int64_t *const aord_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (aord_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_aord_num",
									 "aord_num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.aord_num > 0);
	*aord_num = ctx->jastrow.aord_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_bord_num_device(qmckl_context_device context,
								  int64_t *const bord_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (bord_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_bord_num",
									 "aord_num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 1;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.bord_num > 0);
	*bord_num = ctx->jastrow.bord_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_cord_num_device(qmckl_context_device context,
								  int64_t *const cord_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (cord_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_cord_num",
									 "aord_num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 2;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.cord_num > 0);
	*cord_num = ctx->jastrow.cord_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_type_nucl_num_device(qmckl_context_device context,
									   int64_t *const type_nucl_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (type_nucl_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_type_nucl_num",
									 "type_nucl_num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 3;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.type_nucl_num > 0);
	*type_nucl_num = ctx->jastrow.type_nucl_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_type_nucl_vector_device(qmckl_context_device context,
										  int64_t *const type_nucl_vector,
										  const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (type_nucl_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_type_nucl_vector",
									 "type_nucl_vector is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 4;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.type_nucl_vector != NULL);
	if (size_max < ctx->jastrow.type_nucl_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_type_nucl_vector",
			"Array too small. Expected jastrow.type_nucl_num");
	}

	qmckl_memcpy_D2D(context, type_nucl_vector, ctx->jastrow.type_nucl_vector,
					 ctx->jastrow.type_nucl_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_a_vector_device(qmckl_context_device context,
								  double *const a_vector,
								  const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (a_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_a_vector",
									 "a_vector is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 5;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.a_vector != NULL);
	int64_t sze = (ctx->jastrow.aord_num + 1) * ctx->jastrow.type_nucl_num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE, "qmckl_get_jastrow_a_vector",
			"Array too small. Expected (ctx->jastrow.aord_num + "
			"1)*ctx->jastrow.type_nucl_num");
	}
	qmckl_memcpy_D2D(context, a_vector, ctx->jastrow.a_vector,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_b_vector_device(qmckl_context_device context,
								  double *const b_vector,
								  const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (b_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_b_vector",
									 "b_vector is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 6;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.b_vector != NULL);
	int64_t sze = ctx->jastrow.bord_num + 1;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE, "qmckl_get_jastrow_b_vector",
			"Array too small. Expected (ctx->jastrow.bord_num + 1)");
	}
	qmckl_memcpy_D2D(context, b_vector, ctx->jastrow.b_vector,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_c_vector_device(qmckl_context_device context,
								  double *const c_vector,
								  const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return (char)0;
	}

	if (c_vector == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_c_vector",
									 "c_vector is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 7;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->jastrow.c_vector != NULL);

	int64_t dim_c_vector;
	qmckl_exit_code_device rc =
		qmckl_get_jastrow_dim_c_vector_device(context, &dim_c_vector);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t sze = dim_c_vector * ctx->jastrow.type_nucl_num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_c_vector",
									 "Array too small. Expected dim_c_vector * "
									 "jastrow.type_nucl_num");
	}
	qmckl_memcpy_D2D(context, c_vector, ctx->jastrow.c_vector,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_rescale_factor_ee_device(const qmckl_context_device context,
										   double *const rescale_factor_ee) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (rescale_factor_ee == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_rescale_factor_ee",
									 "rescale_factor_ee is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 8;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}
	assert(ctx->jastrow.rescale_factor_ee > 0.0);

	*rescale_factor_ee = ctx->jastrow.rescale_factor_ee;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_rescale_factor_en_device(const qmckl_context_device context,
										   double *const rescale_factor_en,
										   const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (rescale_factor_en == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_jastrow_rescale_factor_en",
									 "rescale_factor_en is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 9;

	if ((ctx->jastrow.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	if (size_max < ctx->jastrow.type_nucl_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_rescale_factor_en",
									 "Array to small");
	}

	assert(ctx->jastrow.rescale_factor_en != NULL);
	// Originally :
	/*
	for (int64_t i = 0; i < ctx->jastrow.type_nucl_num; ++i) {
		rescale_factor_en[i] = ctx->jastrow.rescale_factor_en[i];
	}
	*/
	qmckl_memcpy_D2D(context, rescale_factor_en, ctx->jastrow.rescale_factor_en,
					 ctx->jastrow.type_nucl_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_dim_c_vector_device(qmckl_context_device context,
									  int64_t *const dim_c_vector) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	*dim_c_vector = ctx->jastrow.dim_c_vector;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_asymp_jasb_device(qmckl_context_device context,
									double *const asymp_jasb,
									const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_asymp_jasb_device",
									 NULL);
	}

	/* Provided in finalize_jastrow */
	/*
	qmckl_exit_code rc;
	rc = qmckl_provide_jastrow_asymp_jasb(context);
	if(rc != QMCKL_SUCCESS) return rc;
	*/

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = 2;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_asymp_jasb_device",
									 "Array too small. Expected 2");
	}
	qmckl_memcpy_D2D(context, asymp_jasb, ctx->jastrow.asymp_jasb,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_asymp_jasa_device(qmckl_context_device context,
									double *const asymp_jasa,
									const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_asymp_jasa_device",
									 NULL);
	}

	/* Provided in finalize_jastrow */
	/*
	qmckl_exit_code rc;
	rc = qmckl_provide_jastrow_asymp_jasa(context);
	if(rc != QMCKL_SUCCESS) return rc;
	*/

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->jastrow.type_nucl_num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_asymp_jasa",
									 "Array too small. Expected nucleus.num");
	}
	qmckl_memcpy_D2D(context, asymp_jasa, ctx->jastrow.asymp_jasa,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// PROVIDE
//**********

// Finalize provides

qmckl_exit_code_device
qmckl_provide_jastrow_asymp_jasa_device(qmckl_context_device context) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_asymp_jasa_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_asymp_jasa_device",
									 NULL);
	}

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.asymp_jasa_date) {

		/* Allocate array */
		if (ctx->jastrow.asymp_jasa == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			double *asymp_jasa = (double *)qmckl_malloc_device(
				context, ctx->jastrow.type_nucl_num * sizeof(double));

			if (asymp_jasa == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "qmckl_asymp_jasa", NULL);
			}
			ctx->jastrow.asymp_jasa = asymp_jasa;
		}

		rc = qmckl_compute_jastrow_asymp_jasa_device(
			context, ctx->jastrow.aord_num, ctx->jastrow.type_nucl_num,
			ctx->jastrow.a_vector, ctx->jastrow.rescale_factor_en,
			ctx->jastrow.asymp_jasa);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.asymp_jasa_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_asymp_jasb_device(qmckl_context_device context) {

	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_asymp_jasb_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_asymp_jasb_device",
									 NULL);
	}

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.asymp_jasb_date) {

		/* Allocate array */
		if (ctx->jastrow.asymp_jasb == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			double *asymp_jasb =
				(double *)qmckl_malloc_device(context, 2 * sizeof(double));

			if (asymp_jasb == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "qmckl_asymp_jasb_device", NULL);
			}
			ctx->jastrow.asymp_jasb = asymp_jasb;
		}

		rc = qmckl_compute_jastrow_asymp_jasb_device(
			context, ctx->jastrow.bord_num, ctx->jastrow.b_vector,
			ctx->jastrow.rescale_factor_ee, ctx->jastrow.asymp_jasb);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.asymp_jasb_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Total Jastrow
qmckl_exit_code_device
qmckl_provide_jastrow_value_device(qmckl_context_device context) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_value", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_value", NULL);
	}

	rc = qmckl_provide_jastrow_factor_ee_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_jastrow_factor_en_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_jastrow_factor_een_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.value_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.value != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.value);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_value",
						"Unable to free ctx->jastrow.value");
				}
				ctx->jastrow.value = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.value == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.walker.num * sizeof(double);
			double *value =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (value == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_value", NULL);
			}
			ctx->jastrow.value = value;
		}

		rc = qmckl_compute_jastrow_value_device(
			context, ctx->electron.walker.num, ctx->jastrow.factor_ee,
			ctx->jastrow.factor_en, ctx->jastrow.factor_een,
			ctx->jastrow.value);

		ctx->jastrow.value_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_gl_device(qmckl_context_device context) {

	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_gl_device", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_gl", NULL);
	}

	rc = qmckl_provide_jastrow_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_jastrow_factor_ee_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_jastrow_factor_en_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_jastrow_factor_een_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.gl_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.gl != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.gl);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_gl",
						"Unable to free ctx->jastrow.gl");
				}
				ctx->jastrow.gl = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.gl == NULL) {
			double *gl = (double *)qmckl_malloc_device(
				context, ctx->electron.walker.num * ctx->electron.num * 4 *
							 sizeof(double));

			if (gl == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_gl_device", NULL);
			}
			ctx->jastrow.gl = gl;
		}

		rc = qmckl_compute_jastrow_gl_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->jastrow.value, ctx->jastrow.factor_ee_deriv_e,
			ctx->jastrow.factor_en_deriv_e, ctx->jastrow.factor_een_deriv_e,
			ctx->jastrow.gl);

		ctx->jastrow.gl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_ee_device(qmckl_context_device context) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_factor_ee", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_factor_ee", NULL);
	}

	rc = qmckl_provide_ee_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_ee_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_ee != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.factor_ee);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_ee",
						"Unable to free ctx->jastrow.factor_ee");
				}
				ctx->jastrow.factor_ee = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_ee == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			double *factor_ee = (double *)qmckl_malloc_device(
				context, ctx->electron.walker.num * sizeof(double));

			if (factor_ee == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_ee", NULL);
			}
			ctx->jastrow.factor_ee = factor_ee;
		}

		rc = qmckl_compute_jastrow_factor_ee_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->electron.up_num, ctx->jastrow.bord_num, ctx->jastrow.b_vector,
			ctx->jastrow.ee_distance_rescaled, ctx->jastrow.asymp_jasb,
			ctx->jastrow.factor_ee);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_ee_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_en_device(qmckl_context_device context) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_factor_en", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_factor_en", NULL);
	}

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_en_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Provided in finalize_jastrow */
	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_en_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_en != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.factor_en);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_en",
						"Unable to free ctx->jastrow.factor_en");
				}
				ctx->jastrow.factor_en = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_en == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.walker.num * sizeof(double);
			double *factor_en =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (factor_en == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_en", NULL);
			}
			ctx->jastrow.factor_en = factor_en;
		}

		rc = qmckl_compute_jastrow_factor_en_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.type_nucl_num,
			ctx->jastrow.type_nucl_vector, ctx->jastrow.aord_num,
			ctx->jastrow.a_vector, ctx->jastrow.en_distance_rescaled,
			ctx->jastrow.asymp_jasa, ctx->jastrow.factor_en);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_en_date = ctx->date;
	}
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_factor_en_deriv_e_device(qmckl_context_device context) {

	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_provide_jastrow_factor_en_deriv_e_device", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(
			context, QMCKL_NOT_PROVIDED_DEVICE,
			"qmckl_provide_jastrow_factor_en_deriv_e_device", NULL);
	}

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_en_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_en_distance_rescaled_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_en_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_en_deriv_e != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.factor_en_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_en_deriv_e",
						"Unable to free ctx->jastrow.factor_en_deriv_e");
				}
				ctx->jastrow.factor_en_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_en_deriv_e == NULL) {

			double *factor_en_deriv_e = (double *)qmckl_malloc_device(
				context, ctx->electron.walker.num * 4 * ctx->electron.num *
							 sizeof(double));

			if (factor_en_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_en_deriv_e_device", NULL);
			}
			ctx->jastrow.factor_en_deriv_e = factor_en_deriv_e;
		}

		rc = qmckl_compute_jastrow_factor_en_deriv_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.type_nucl_num,
			ctx->jastrow.type_nucl_vector, ctx->jastrow.aord_num,
			ctx->jastrow.a_vector, ctx->jastrow.en_distance_rescaled,
			ctx->jastrow.en_distance_rescaled_deriv_e,
			ctx->jastrow.factor_en_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_en_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_een_device(qmckl_context_device context) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_een_rescaled_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_jastrow_c_vector_full_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_lkpm_combined_index_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if tmp_c is provided */
	rc = qmckl_provide_tmp_c_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_een_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_een != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.factor_een);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_een",
						"Unable to free ctx->jastrow.factor_een");
				}
				ctx->jastrow.factor_een = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_een == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.walker.num * sizeof(double);
			double *factor_een =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (factor_een == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_een", NULL);
			}
			ctx->jastrow.factor_een = factor_een;
		}

		rc = qmckl_compute_jastrow_factor_een_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.cord_num, ctx->jastrow.dim_c_vector,
			ctx->jastrow.c_vector_full, ctx->jastrow.lkpm_combined_index,
			ctx->jastrow.tmp_c, ctx->jastrow.een_rescaled_n,
			ctx->jastrow.factor_een);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_een_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_factor_een_deriv_e_device(qmckl_context_device context) {

	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_een_rescaled_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_een_rescaled_e_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_een_rescaled_n_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_jastrow_c_vector_full_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_lkpm_combined_index_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if tmp_c is provided */
	rc = qmckl_provide_tmp_c_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if dtmp_c is provided */
	rc = qmckl_provide_dtmp_c_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_een_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_een_deriv_e != NULL) {
				rc =
					qmckl_free_device(context, ctx->jastrow.factor_een_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_een_deriv_e",
						"Unable to free ctx->jastrow.factor_een_deriv_e");
				}
				ctx->jastrow.factor_een_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_een_deriv_e == NULL) {

			double *factor_een_deriv_e = (double *)qmckl_malloc_device(
				context, 4 * ctx->electron.num * ctx->electron.walker.num *
							 sizeof(double));

			if (factor_een_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_een_deriv_e_device", NULL);
			}
			ctx->jastrow.factor_een_deriv_e = factor_een_deriv_e;
		}

		rc = qmckl_compute_jastrow_factor_een_deriv_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.cord_num, ctx->jastrow.dim_c_vector,
			ctx->jastrow.c_vector_full, ctx->jastrow.lkpm_combined_index,
			ctx->jastrow.tmp_c, ctx->jastrow.dtmp_c,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.een_rescaled_n_deriv_e,
			ctx->jastrow.factor_een_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_een_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_dtmp_c_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;
	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	rc = qmckl_provide_een_rescaled_e_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.dtmp_c_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.dtmp_c != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.dtmp_c);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_dtmp_c_device",
						"Unable to free ctx->jastrow.dtmp_c");
				}
				ctx->jastrow.dtmp_c = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.dtmp_c == NULL) {

			double *dtmp_c = (double *)qmckl_malloc_device(
				context, (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1) *
							 4 * ctx->electron.num * ctx->nucleus.num *
							 ctx->electron.walker.num * sizeof(double));

			if (dtmp_c == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_dtmp_c_device", NULL);
			}
			ctx->jastrow.dtmp_c = dtmp_c;
		}

		rc = qmckl_compute_dtmp_c_device(
			context, ctx->jastrow.cord_num, ctx->electron.num, ctx->nucleus.num,
			ctx->electron.walker.num, ctx->jastrow.een_rescaled_e_deriv_e,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.dtmp_c);

		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.dtmp_c_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Distances
qmckl_exit_code_device
qmckl_provide_ee_distance_rescaled_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Compute if necessary */
	if (ctx->electron.walker.point.date >
		ctx->jastrow.ee_distance_rescaled_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.ee_distance_rescaled != NULL) {
				qmckl_exit_code_device rc = qmckl_free_device(
					context, ctx->jastrow.ee_distance_rescaled);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_ee_distance_rescaled",
						"Unable to free "
						"ctx->jastrow.ee_distance_rescaled");
				}
				ctx->jastrow.ee_distance_rescaled = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.ee_distance_rescaled == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.num * ctx->electron.num *
							ctx->electron.walker.num * sizeof(double);
			double *ee_distance_rescaled =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (ee_distance_rescaled == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_ee_distance_rescaled", NULL);
			}
			ctx->jastrow.ee_distance_rescaled = ee_distance_rescaled;
		}

		qmckl_exit_code_device rc = qmckl_compute_ee_distance_rescaled_device(
			context, ctx->electron.num, ctx->jastrow.rescale_factor_ee,
			ctx->electron.walker.num, ctx->electron.walker.point.coord.data,
			ctx->jastrow.ee_distance_rescaled);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.ee_distance_rescaled_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_en_distance_rescaled_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!(ctx->nucleus.provided)) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	/* Compute if necessary */
	if (ctx->electron.walker.point.date >
		ctx->jastrow.en_distance_rescaled_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.en_distance_rescaled != NULL) {
				qmckl_exit_code_device rc = qmckl_free_device(
					context, ctx->jastrow.en_distance_rescaled);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_en_distance_rescaled",
						"Unable to free "
						"ctx->jastrow.en_distance_rescaled");
				}
				ctx->jastrow.en_distance_rescaled = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.en_distance_rescaled == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.num * ctx->nucleus.num *
							ctx->electron.walker.num * sizeof(double);
			double *en_distance_rescaled =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (en_distance_rescaled == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_en_distance_rescaled", NULL);
			}
			ctx->jastrow.en_distance_rescaled = en_distance_rescaled;
		}

		qmckl_exit_code_device rc = qmckl_compute_en_distance_rescaled_device(
			context, ctx->electron.num, ctx->nucleus.num,
			ctx->jastrow.type_nucl_num, ctx->jastrow.type_nucl_vector,
			ctx->jastrow.rescale_factor_en, ctx->electron.walker.num,
			ctx->electron.walker.point.coord.data, ctx->nucleus.coord.data,
			ctx->jastrow.en_distance_rescaled);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.en_distance_rescaled_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_provide_en_distance_rescaled_deriv_e_device(
	qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->nucleus.provided) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	/* Compute if necessary */
	if (ctx->electron.walker.point.date >
		ctx->jastrow.en_distance_rescaled_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.en_distance_rescaled_deriv_e != NULL) {
				qmckl_exit_code_device rc = qmckl_free_device(
					context, ctx->jastrow.en_distance_rescaled_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc,
						"qmckl_provide_en_distance_rescaled_deriv_e",
						"Unable to free "
						"ctx->jastrow.en_distance_rescaled_deriv_e");
				}
				ctx->jastrow.en_distance_rescaled_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.en_distance_rescaled_deriv_e == NULL) {

			double *en_distance_rescaled_deriv_e =
				(double *)qmckl_malloc_device(
					context, 4 * ctx->electron.num * ctx->nucleus.num *
								 ctx->electron.walker.num * sizeof(double));

			if (en_distance_rescaled_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_en_distance_rescaled_deriv_e_device", NULL);
			}
			ctx->jastrow.en_distance_rescaled_deriv_e =
				en_distance_rescaled_deriv_e;
		}

		qmckl_exit_code_device rc =
			qmckl_compute_en_distance_rescaled_deriv_e_device(
				context, ctx->electron.num, ctx->nucleus.num,
				ctx->jastrow.type_nucl_num, ctx->jastrow.type_nucl_vector,
				ctx->jastrow.rescale_factor_en, ctx->electron.walker.num,
				ctx->electron.walker.point.coord.data, ctx->nucleus.coord.data,
				ctx->jastrow.en_distance_rescaled_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.en_distance_rescaled_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_een_rescaled_e_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if ee distance is provided */
	qmckl_exit_code_device rc = qmckl_provide_ee_distance_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.een_rescaled_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.een_rescaled_e != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.een_rescaled_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_een_rescaled_e",
						"Unable to free ctx->jastrow.een_rescaled_e");
				}
				ctx->jastrow.een_rescaled_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.een_rescaled_e == NULL) {
			double *een_rescaled_e = (double *)qmckl_malloc_device(
				context, ctx->electron.num * ctx->electron.num *
							 ctx->electron.walker.num *
							 (ctx->jastrow.cord_num + 1) * sizeof(double));

			if (een_rescaled_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_een_rescaled_e", NULL);
			}
			ctx->jastrow.een_rescaled_e = een_rescaled_e;
		}

		rc = qmckl_compute_een_rescaled_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->jastrow.cord_num, ctx->jastrow.rescale_factor_ee,
			ctx->electron.ee_distance, ctx->jastrow.een_rescaled_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.een_rescaled_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_een_rescaled_n_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if ee distance is provided */
	qmckl_exit_code_device rc = qmckl_provide_en_distance_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.een_rescaled_n_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.een_rescaled_n != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.een_rescaled_n);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_een_rescaled_n",
						"Unable to free ctx->jastrow.een_rescaled_n");
				}
				ctx->jastrow.een_rescaled_n = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.een_rescaled_n == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.num * ctx->nucleus.num *
							ctx->electron.walker.num *
							(ctx->jastrow.cord_num + 1) * sizeof(double);
			double *een_rescaled_n =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (een_rescaled_n == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_een_rescaled_n", NULL);
			}
			ctx->jastrow.een_rescaled_n = een_rescaled_n;
		}

		rc = qmckl_compute_een_rescaled_n_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.type_nucl_num,
			ctx->jastrow.type_nucl_vector, ctx->jastrow.cord_num,
			ctx->jastrow.rescale_factor_en, ctx->electron.en_distance,
			ctx->jastrow.een_rescaled_n);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.een_rescaled_n_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_c_vector_full_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.c_vector_full_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.c_vector_full != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.c_vector_full);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_c_vector_full",
						"Unable to free "
						"ctx->jastrow.c_vector_full");
				}
				ctx->jastrow.c_vector_full = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.c_vector_full == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size =
				ctx->jastrow.dim_c_vector * ctx->nucleus.num * sizeof(double);
			double *c_vector_full =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (c_vector_full == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_c_vector_full", NULL);
			}
			ctx->jastrow.c_vector_full = c_vector_full;
		}

		rc = qmckl_compute_c_vector_full_device(
			context, ctx->nucleus.num, ctx->jastrow.dim_c_vector,
			ctx->jastrow.type_nucl_num, ctx->jastrow.type_nucl_vector,
			ctx->jastrow.c_vector, ctx->jastrow.c_vector_full);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.c_vector_full_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_lkpm_combined_index_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.lkpm_combined_index_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.lkpm_combined_index != NULL) {
				rc = qmckl_free_device(context,
									   ctx->jastrow.lkpm_combined_index);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_ee",
						"Unable to free "
						"ctx->jastrow.lkpm_combined_index");
				}
				ctx->jastrow.lkpm_combined_index = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.lkpm_combined_index == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = 4 * ctx->jastrow.dim_c_vector * sizeof(int64_t);
			int64_t *lkpm_combined_index =
				(int64_t *)qmckl_malloc_device(context, mem_info.size);

			if (lkpm_combined_index == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_lkpm_combined_index", NULL);
			}
			ctx->jastrow.lkpm_combined_index = lkpm_combined_index;
		}

		rc = qmckl_compute_lkpm_combined_index_device(
			context, ctx->jastrow.cord_num, ctx->jastrow.dim_c_vector,
			ctx->jastrow.lkpm_combined_index);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.lkpm_combined_index_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_tmp_c_device(qmckl_context_device context) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	rc = qmckl_provide_een_rescaled_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.tmp_c_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.tmp_c != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.tmp_c);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_tmp_c",
						"Unable to free ctx->jastrow.tmp_c");
				}
				ctx->jastrow.tmp_c = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.tmp_c == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = (ctx->jastrow.cord_num) *
							(ctx->jastrow.cord_num + 1) * ctx->electron.num *
							ctx->nucleus.num * ctx->electron.walker.num *
							sizeof(double);
			double *tmp_c =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (tmp_c == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "qmckl_provide_tmp_c", NULL);
			}
			ctx->jastrow.tmp_c = tmp_c;
		}

		rc = qmckl_compute_tmp_c_device(
			context, ctx->jastrow.cord_num, ctx->electron.num, ctx->nucleus.num,
			ctx->electron.walker.num, ctx->jastrow.een_rescaled_e,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.tmp_c);

		ctx->jastrow.tmp_c_date = ctx->date;
	}
	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron/nucleus deriv
qmckl_exit_code_device
qmckl_provide_een_rescaled_e_deriv_e_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if ee distance is provided */
	qmckl_exit_code_device rc = qmckl_provide_een_rescaled_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.een_rescaled_e_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.een_rescaled_e_deriv_e != NULL) {
				rc = qmckl_free_device(context,
									   ctx->jastrow.een_rescaled_e_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_een_rescaled_e_deriv_e",
						"Unable to free ctx->jastrow.een_rescaled_e_deriv_e");
				}
				ctx->jastrow.een_rescaled_e_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.een_rescaled_e_deriv_e == NULL) {

			qmckl_memory_info_struct_device mem_info =
				qmckl_memory_info_struct_zero_device;
			mem_info.size = ctx->electron.num * 4 * ctx->electron.num *
							ctx->electron.walker.num *
							(ctx->jastrow.cord_num + 1) * sizeof(double);
			double *een_rescaled_e_deriv_e =
				(double *)qmckl_malloc_device(context, mem_info.size);

			if (een_rescaled_e_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_een_rescaled_e_deriv_e", NULL);
			}
			ctx->jastrow.een_rescaled_e_deriv_e = een_rescaled_e_deriv_e;
		}

		rc = qmckl_compute_jastrow_factor_een_rescaled_e_deriv_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->jastrow.cord_num, ctx->jastrow.rescale_factor_ee,
			ctx->electron.walker.point.coord.data, ctx->electron.ee_distance,
			ctx->jastrow.een_rescaled_e, ctx->jastrow.een_rescaled_e_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.een_rescaled_e_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_een_rescaled_n_deriv_e_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Check if ee distance is provided */
	qmckl_exit_code_device rc = qmckl_provide_en_distance_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if ee distance is provided */
	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.een_rescaled_n_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.een_rescaled_n_deriv_e != NULL) {
				rc = qmckl_free_device(context,
									   ctx->jastrow.een_rescaled_n_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc,
						"qmckl_provide_een_rescaled_n_deriv_e_device",
						"Unable to free ctx->jastrow.een_rescaled_n_deriv_e");
				}
				ctx->jastrow.een_rescaled_n_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.een_rescaled_n_deriv_e == NULL) {

			double *een_rescaled_n_deriv_e = (double *)qmckl_malloc_device(
				context, ctx->electron.num * 4 * ctx->nucleus.num *
							 ctx->electron.walker.num *
							 (ctx->jastrow.cord_num + 1) * sizeof(double));

			if (een_rescaled_n_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_een_rescaled_n_deriv_e_device", NULL);
			}
			ctx->jastrow.een_rescaled_n_deriv_e = een_rescaled_n_deriv_e;
		}

		rc = qmckl_compute_jastrow_factor_een_rescaled_n_deriv_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.type_nucl_num,
			ctx->jastrow.type_nucl_vector, ctx->jastrow.cord_num,
			ctx->jastrow.rescale_factor_en,
			ctx->electron.walker.point.coord.data, ctx->nucleus.coord.data,
			ctx->electron.en_distance, ctx->jastrow.een_rescaled_n,
			ctx->jastrow.een_rescaled_n_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.een_rescaled_n_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_provide_jastrow_factor_ee_deriv_e_device(qmckl_context_device context) {

	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_jastrow_factor_ee_deriv_e",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->jastrow.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_jastrow_factor_ee_deriv_e",
									 NULL);
	}

	/* Check if ee rescaled distance is provided */
	rc = qmckl_provide_ee_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Check if ee rescaled distance deriv e is provided */
	rc = qmckl_provide_ee_distance_rescaled_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_ee_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.factor_ee_deriv_e != NULL) {
				rc = qmckl_free_device(context, ctx->jastrow.factor_ee_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc, "qmckl_provide_jastrow_factor_ee_deriv_e",
						"Unable to free ctx->jastrow.factor_ee_deriv_e");
				}
				ctx->jastrow.factor_ee_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.factor_ee_deriv_e == NULL) {

			double *factor_ee_deriv_e = (double *)qmckl_malloc_device(
				context, ctx->electron.walker.num * 4 * ctx->electron.num *
							 sizeof(double));

			if (factor_ee_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_jastrow_factor_ee_deriv_e", NULL);
			}
			ctx->jastrow.factor_ee_deriv_e = factor_ee_deriv_e;
		}

		rc = qmckl_compute_jastrow_factor_ee_deriv_e_device(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->electron.up_num, ctx->jastrow.bord_num, ctx->jastrow.b_vector,
			ctx->jastrow.ee_distance_rescaled,
			ctx->jastrow.ee_distance_rescaled_deriv_e,
			ctx->jastrow.factor_ee_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.factor_ee_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_provide_ee_distance_rescaled_deriv_e_device(
	qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	/* Compute if necessary */
	if (ctx->electron.walker.point.date >
		ctx->jastrow.ee_distance_rescaled_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			if (ctx->jastrow.ee_distance_rescaled_deriv_e != NULL) {
				qmckl_exit_code_device rc = qmckl_free_device(
					context, ctx->jastrow.ee_distance_rescaled_deriv_e);
				if (rc != QMCKL_SUCCESS_DEVICE) {
					return qmckl_failwith_device(
						context, rc,
						"qmckl_provide_ee_distance_rescaled_deriv_e_device",
						"Unable to free "
						"ctx->jastrow.ee_distance_rescaled_deriv_e");
				}
				ctx->jastrow.ee_distance_rescaled_deriv_e = NULL;
			}
		}

		/* Allocate array */
		if (ctx->jastrow.ee_distance_rescaled_deriv_e == NULL) {
			double *ee_distance_rescaled_deriv_e =
				(double *)qmckl_malloc_device(
					context, 4 * ctx->electron.num * ctx->electron.num *
								 ctx->electron.walker.num * sizeof(double));

			if (ee_distance_rescaled_deriv_e == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_ee_distance_rescaled_deriv_e_device", NULL);
			}
			ctx->jastrow.ee_distance_rescaled_deriv_e =
				ee_distance_rescaled_deriv_e;
		}

		qmckl_exit_code_device rc =
			qmckl_compute_ee_distance_rescaled_deriv_e_device(
				context, ctx->electron.num, ctx->jastrow.rescale_factor_ee,
				ctx->electron.walker.num, ctx->electron.walker.point.coord.data,
				ctx->jastrow.ee_distance_rescaled_deriv_e);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->jastrow.ee_distance_rescaled_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GETTERS (for computes)
//**********

// Total Jastrow
qmckl_exit_code_device
qmckl_get_jastrow_value_device(qmckl_context_device context,
							   double *const value, const int64_t size_max) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_value", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	rc = qmckl_provide_jastrow_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t sze = ctx->electron.walker.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_value",
									 "Array too small. Expected walker.num");
	}
	qmckl_memcpy_D2D(context, value, ctx->jastrow.value, sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_jastrow_gl_device(qmckl_context_device context,
												   double *const gl,
												   const int64_t size_max) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_champ_gl_device", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	rc = qmckl_provide_jastrow_gl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t sze = 4 * ctx->electron.walker.num * ctx->electron.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_champ_gl_device",
			"Array too small. Expected walker.num * electron.num * 4");
	}
	qmckl_memcpy_D2D(context, gl, ctx->jastrow.gl, sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron component
qmckl_exit_code_device
qmckl_get_jastrow_factor_ee_device(qmckl_context_device context,
								   double *const factor_ee,
								   const int64_t size_max) {
	qmckl_exit_code_device rc;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_factor_ee", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	rc = qmckl_provide_jastrow_factor_ee_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t sze = ctx->electron.walker.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_factor_ee",
									 "Array too small. Expected walker.num");
	}
	qmckl_memcpy_D2D(context, factor_ee, ctx->jastrow.factor_ee,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

// Electron/nucleus component
qmckl_exit_code_device
qmckl_get_jastrow_factor_en_device(qmckl_context_device context,
								   double *const factor_en,
								   const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_jastrow_factor_en", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_factor_en_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t sze = ctx->electron.walker.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_factor_en",
									 "Array too small. Expected walker.num");
	}
	qmckl_memcpy_D2D(context, factor_en, ctx->jastrow.factor_en,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_factor_en_deriv_e_device(qmckl_context_device context,
										   double *const factor_en_deriv_e,
										   const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_factor_en_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_en_deriv_e",
			"Array too small. Expected 4*walker.num*elec_num");
	}
	qmckl_memcpy_D2D(context, factor_en_deriv_e, ctx->jastrow.factor_en_deriv_e,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_get_jastrow_factor_een_device(qmckl_context_device context,
									double *const factor_een,
									const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_factor_een_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.walker.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_jastrow_factor_een",
									 "Array too small. Expected walk_num");
	}
	qmckl_memcpy_D2D(context, factor_een, ctx->jastrow.factor_een,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_factor_een_deriv_e_device(qmckl_context_device context,
											double *const factor_een_deriv_e,
											const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_factor_een_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_een_deriv_e_device",
			"Array too small. Expected 4*walk_num*elec_num");
	}
	qmckl_memcpy_D2D(context, factor_een_deriv_e,
					 ctx->jastrow.factor_een_deriv_e, sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_e_deriv_e_device(qmckl_context_device context,
												double *const distance_rescaled,
												const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_een_rescaled_e_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.num * 4 * ctx->electron.num *
				  ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_een_deriv_e",
			"Array too small. Expected ctx->electron.num * 4 * "
			"ctx->electron.num * ctx->electron.walker.num * "
			"(ctx->jastrow.cord_num + 1)");
	}
	qmckl_memcpy_D2D(context, distance_rescaled,
					 ctx->jastrow.een_rescaled_e_deriv_e, sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_n_deriv_e_device(qmckl_context_device context,
												double *const distance_rescaled,
												const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_een_rescaled_n_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.num * 4 * ctx->nucleus.num *
				  ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_een_deriv_e",
			"Array too small. Expected ctx->electron.num * 4 * "
			"ctx->nucleus.num * ctx->electron.walker.num * "
			"(ctx->jastrow.cord_num + 1)");
	}
	qmckl_memcpy_D2D(context, distance_rescaled,
					 ctx->jastrow.een_rescaled_n_deriv_e, sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_dtmp_c_device(qmckl_context_device context,
								double *const dtmp_c) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_c_vector_full_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_dtmp_c_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1) * 4 *
				 ctx->electron.num * ctx->nucleus.num *
				 ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, dtmp_c, ctx->jastrow.dtmp_c,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

// Distances
qmckl_exit_code_device
qmckl_get_jastrow_ee_distance_rescaled_device(qmckl_context_device context,
											  double *const distance_rescaled) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ee_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze =
		ctx->electron.num * ctx->electron.num * ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, distance_rescaled,
					 ctx->jastrow.ee_distance_rescaled, sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device(
	qmckl_context_device context, double *const distance_rescaled_deriv_e) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ee_distance_rescaled_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze =
		4 * ctx->electron.num * ctx->electron.num * ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, distance_rescaled_deriv_e,
					 ctx->jastrow.ee_distance_rescaled_deriv_e,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_electron_en_distance_rescaled_device(qmckl_context_device context,
											   double *distance_rescaled) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_en_distance_rescaled_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze =
		ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, distance_rescaled,
					 ctx->jastrow.en_distance_rescaled, sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_electron_en_distance_rescaled_deriv_e_device(
	qmckl_context_device context, double *distance_rescaled_deriv_e) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_en_distance_rescaled_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze =
		4 * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, distance_rescaled_deriv_e,
					 ctx->jastrow.en_distance_rescaled_deriv_e,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_e_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_een_rescaled_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.num * ctx->electron.num *
				  ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);

	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_een_deriv_e",
			"Array too small. Expected ctx->electron.num * 4 * "
			"ctx->electron.num * ctx->electron.walker.num * "
			"(ctx->jastrow.cord_num + 1)");
	}

	qmckl_memcpy_D2D(context, distance_rescaled, ctx->jastrow.een_rescaled_e,
					 sze * sizeof(double));
#pragma acc kernels deviceptr(distance_rescaled)
	{}
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_n_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_een_rescaled_n_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.num * ctx->nucleus.num *
				  ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_een_deriv_e",
			"Array too small. Expected ctx->electron.num * "
			"ctx->nucleus.num * ctx->electron.walker.num * "
			"(ctx->jastrow.cord_num + 1)");
	}
	qmckl_memcpy_D2D(context, distance_rescaled, ctx->jastrow.een_rescaled_n,
					 sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_tmp_c_device(qmckl_context_device context,
							   double *const tmp_c) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_c_vector_full_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	rc = qmckl_provide_tmp_c_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1) *
				 ctx->electron.num * ctx->nucleus.num *
				 ctx->electron.walker.num;
	qmckl_memcpy_D2D(context, tmp_c, ctx->jastrow.tmp_c, sze * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

// Misc

// This simply computes a scalar value by side effect
qmckl_exit_code_device
qmckl_compute_dim_c_vector_device(const qmckl_context_device context,
								  const int64_t cord_num,
								  int64_t *const dim_c_vector) {

	int lmax;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (cord_num < 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	*dim_c_vector = 0;

	for (int p = 2; p <= cord_num; ++p) {
		for (int k = p - 1; k >= 0; --k) {
			if (k != 0) {
				lmax = p - k;
			} else {
				lmax = p - k - 2;
			}
			for (int l = lmax; l >= 0; --l) {
				if (((p - k - l) & 1) == 1)
					continue;
				*dim_c_vector = *dim_c_vector + 1;
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_jastrow_factor_ee_deriv_e_device(qmckl_context_device context,
										   double *const factor_ee_deriv_e,
										   const int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_jastrow_factor_ee_deriv_e_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_jastrow_factor_ee_deriv_e_device",
			"Array too small. Expected 4*walk_num*elec_num");
	}

	qmckl_memcpy_D2D(context, factor_ee_deriv_e, ctx->jastrow.factor_ee_deriv_e,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}
