#include "../include/qmckl_point.h"

qmckl_exit_code_device qmckl_set_point_device(qmckl_context_device context,
											  char transp, int64_t num,
											  double *coord, int64_t size_max) {

	size_t device_id = qmckl_get_device_id(context);
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	if (size_max < 3 * num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_set_point_device",
									 "Array too small");
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_point_device",
									 "transp should be 'N' or 'T'");
	}

	if (coord == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_point_device",
									 "coord is a NULL pointer");
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc;
	if (num != ctx->point.num) {

		if (ctx->point.coord.data != NULL) {
			rc = qmckl_matrix_free_device(context, &(ctx->point.coord));
			assert(rc == QMCKL_SUCCESS_DEVICE);
		}

		ctx->point.coord = qmckl_matrix_alloc_device(context, num, 3);
		if (ctx->point.coord.data == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
										 "qmckl_set_point", NULL);
		}
	};

	ctx->point.num = num;

	double *a = ctx->point.coord.data;
	int size_0 = ctx->point.coord.size[0];
	if (transp == 'T') {
#pragma acc kernels deviceptr(a, coord)
#pragma omp target is_device_ptr(a, coord)
{
		for (int64_t i = 0; i < 3 * num; ++i) {
			a[i] = coord[i];
		}
}
	} else {
#pragma acc kernels deviceptr(a, coord)
#pragma omp target is_device_ptr(a, coord)
{
		for (int64_t i = 0; i < num; ++i) {
			a[i] = coord[3 * i];
			a[i + size_0] = coord[3 * i + 1];
			a[i + 2 * size_0] = coord[3 * i + 2];
		}
}
	}

	/* Increment the date of the context */
	rc = qmckl_context_touch_device(context);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	return QMCKL_SUCCESS_DEVICE;
}
