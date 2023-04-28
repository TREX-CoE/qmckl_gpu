#include "../include/qmckl_nucleus.h"

/* Provided check  */

bool qmckl_nucleus_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->nucleus.provided;
}

//**********
// SETTERS
//**********

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num) {
	int32_t mask = 1 << 0;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	if (num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_nucleus_num", "num <= 0");
	}

	ctx->nucleus.num = num;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_nucleus_coord_device(qmckl_context_device context, char transp,
							   double *coord, int64_t size_max) {
	int32_t mask = 1 << 2;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_exit_code_device rc;

	const int64_t nucl_num = (int64_t)ctx->nucleus.num;

	if (ctx->nucleus.coord.data != NULL) {
		rc = qmckl_matrix_free_device(context, &(ctx->nucleus.coord));
		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	ctx->nucleus.coord = qmckl_matrix_alloc_device(context, nucl_num, 3);

	if (ctx->nucleus.coord.data == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_nucleus_coord", NULL);
	}

	if (size_max < 3 * nucl_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_set_nucleus_coord",
									 "Array too small");
	}

	if (transp == 'N') {
		qmckl_matrix_device At;
		At = qmckl_matrix_alloc_device(context, 3, nucl_num);
		rc = qmckl_matrix_of_double_device(context, coord, 3 * nucl_num, &At);
		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
		rc = qmckl_transpose_device(context, At, ctx->nucleus.coord);
	} else {
		rc = qmckl_matrix_of_double_device(context, coord, nucl_num * 3,
										   &(ctx->nucleus.coord));
	}
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

/* Sets the nuclear charges of all the atoms. */
qmckl_exit_code_device
qmckl_set_nucleus_charge_device(qmckl_context_device context, double *charge,
								int64_t size_max) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_charge_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (charge == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "charge is a null pointer");
	}

	int64_t num;
	qmckl_exit_code_device rc;

	rc = qmckl_get_nucleus_num_device(context, &num);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	if (num > size_max) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "Array too small");
	}

	ctx->nucleus.charge = qmckl_vector_alloc_device(context, num);
	rc = qmckl_vector_of_double_device(context, charge, num,
									   &(ctx->nucleus.charge));

	if (rc != QMCKL_SUCCESS_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "Error in vector->double* conversion");
	}

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// SETTERS
//**********

qmckl_exit_code_device
qmckl_get_nucleus_num_device(const qmckl_context_device context,
							 int64_t *const num) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_nucleus_num",
									 "num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->nucleus.uninitialized & mask) != 0) {
		*num = (int64_t)0;
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_nucleus_num",
									 "nucleus data is not provided");
	}

	assert(ctx->nucleus.num >= (int64_t)0);
	*num = ctx->nucleus.num;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_nucleus_coord_device(const qmckl_context_device context,
							   const char transp, double *const coord,
							   const int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_nucleus_coord_device",
									 "transp should be 'N' or 'T'");
	}

	if (coord == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_nucleus_coord_device",
									 "coord is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 2;

	if ((ctx->nucleus.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_nucleus_coord_device",
									 "nucleus data is not provided");
	}

	assert(ctx->nucleus.coord.data != NULL);

	qmckl_exit_code_device rc;

	if (transp == 'N') {
		qmckl_matrix_device At =
			qmckl_matrix_alloc_device(context, 3, ctx->nucleus.coord.size[0]);
		rc = qmckl_transpose_device(context, ctx->nucleus.coord, At);
		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;

		// Copy content of At onto coord
		// rc = qmckl_double_of_matrix_device(context, At, coord, size_max);
		qmckl_memcpy_D2D(context, coord, At.data,
						 At.size[0] * At.size[1] * sizeof(double));
		qmckl_matrix_free_device(context, &At);
	} else {
		// Copy content of ctx->nucleus.coord onto coord
		// rc = qmckl_double_of_matrix_device(context, ctx->nucleus.coord,
		// coord, 							size_max);
		qmckl_memcpy_D2D(context, coord, ctx->nucleus.coord.data,
						 ctx->nucleus.coord.size[0] *
							 ctx->nucleus.coord.size[1] * sizeof(double));
	}

	return rc;
}

qmckl_exit_code_device
qmckl_get_nucleus_charge_device(const qmckl_context_device context,
								double *const charge, const int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (charge == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_nucleus_charge",
									 "charge is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 1;

	if ((ctx->nucleus.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_nucleus_charge_device",
									 "nucleus data is not provided");
	}

	assert(ctx->nucleus.charge.data != NULL);

	if (ctx->nucleus.num > size_max) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_nucleus_charge",
									 "Array too small");
	}

	// Copy content of ctx->nucleus.charge onto charge
	// rc = qmckl_double_of_vector_device(context, ctx->nucleus.charge, charge,
	// size_max);
	qmckl_memcpy_D2D(context, charge, ctx->nucleus.charge.data,
					 ctx->nucleus.charge.size * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}
