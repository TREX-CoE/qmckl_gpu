#include "../include/qmckl_electron.h"

/* Provided check  */

bool qmckl_electron_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->electron.provided;
}

//**********
// SETTERS
//**********

qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context, int64_t up_num,
							  int64_t down_num) {

	int32_t mask = 1 << 0;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_electron_*", NULL);
	}

	if (up_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_electron_num_device",
									 "up_num <= 0");
	}

	if (down_num < 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_electron_num_device",
									 "down_num < 0");
	}

	ctx->electron.up_num = up_num;
	ctx->electron.down_num = down_num;
	ctx->electron.num = up_num + down_num;

	ctx->electron.uninitialized &= ~mask;
	ctx->electron.provided = (ctx->electron.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_electron_coord_device(qmckl_context_device context, char transp,
								int64_t walk_num, double *coord,
								int64_t size_max) {

	size_t device_id = qmckl_get_device_id(context);
	int32_t mask = 0; // coord can be changed

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int64_t elec_num = ctx->electron.num;

	if (elec_num == 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_electron_coord_device",
									 "elec_num is not set");
	}

	/* Swap pointers */
	qmckl_walker_device tmp = ctx->electron.walker_old;
	ctx->electron.walker_old = ctx->electron.walker;
	ctx->electron.walker = tmp;

	memcpy(&(ctx->point), &(ctx->electron.walker.point),
		   sizeof(qmckl_point_struct_device));

	qmckl_exit_code_device rc;
	rc = qmckl_set_point_device(context, transp, walk_num * elec_num, coord,
								size_max);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->electron.walker.num = walk_num;
	memcpy(&(ctx->electron.walker.point), &(ctx->point),
		   sizeof(qmckl_point_struct_device));

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GETTERS
//**********

qmckl_exit_code_device
qmckl_get_electron_num_device(const qmckl_context_device context,
							  int64_t *const num) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_electron_num",
									 "num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->electron.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED_DEVICE;
	}

	assert(ctx->electron.num > (int64_t)0);
	*num = ctx->electron.num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_electron_coord_device(const qmckl_context_device context,
								const char transp, double *const coord,
								const int64_t size_max) {
	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_electron_coord_device",
									 "transp should be 'N' or 'T'");
	}

	if (coord == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_electron_coord_device",
									 "coord is a null pointer");
	}

	if (size_max <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_get_electron_coord_device",
									 "size_max should be > 0");
	}

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->electron.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_electron_coord_device", NULL);
	}

	assert(ctx->point.num == ctx->electron.walker.point.num);
	assert(ctx->point.coord.data == ctx->electron.walker.point.coord.data);

	return qmckl_get_point_device(context, transp, coord, size_max);
}
