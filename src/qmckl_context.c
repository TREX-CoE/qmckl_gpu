// This file contains qmckl_context_device related functions

#include "../include/qmckl_context.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code_device
qmckl_context_touch_device(const qmckl_context_device context) {
	return qmckl_context_touch_device(context);
}

qmckl_context_device
qmckl_context_check_device(const qmckl_context_device context) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_NULL_CONTEXT_DEVICE;

	const qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	/* Try to access memory */
	if (ctx->tag != VALID_TAG_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	return context;
}

void qmckl_lock_device(qmckl_context_device context) {
	if (context == QMCKL_NULL_CONTEXT_DEVICE)
		return;
	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	errno = 0;
	int rc = pthread_mutex_lock(&(ctx->mutex));
	if (rc != 0) {
		fprintf(stderr, "DEBUG qmckl_lock:%s\n", strerror(rc));
		fflush(stderr);
	}
	assert(rc == 0);
	ctx->lock_count += 1;
}

void qmckl_unlock_device(const qmckl_context_device context) {
	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	int rc = pthread_mutex_unlock(&(ctx->mutex));
	if (rc != 0) {
		fprintf(stderr, "DEBUG qmckl_unlock:%s\n", strerror(rc));
		fflush(stderr);
	}
	assert(rc == 0);
	ctx->lock_count -= 1;
}

//**********
// CONTEXT CREATE
//**********

qmckl_exit_code_device
qmckl_init_electron_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->electron.uninitialized = (1 << 1) - 1;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_init_nucleus_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->nucleus.uninitialized = (1 << 3) - 1;

	/* Default values */

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_init_point_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	memset(&(ctx->point), 0, sizeof(qmckl_point_struct_device));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_init_mo_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->mo_basis.uninitialized = (1 << 2) - 1;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_init_determinant_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->det.uninitialized = (1 << 5) - 1;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->jastrow.uninitialized = (1 << 10) - 1;

	/* Default values */
	ctx->jastrow.aord_num = -1;
	ctx->jastrow.bord_num = -1;
	ctx->jastrow.cord_num = -1;
	ctx->jastrow.type_nucl_num = -1;
	ctx->jastrow.dim_c_vector = -1;

	return QMCKL_SUCCESS_DEVICE;
}


qmckl_context_device qmckl_context_create_device(int device_id) {

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)malloc(
			sizeof(qmckl_context_struct_device));

	if (ctx == NULL) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	/* Set all pointers and values to NULL */
	{ memset(ctx, 0, sizeof(qmckl_context_struct_device)); }

	/* Initialize lock */
	{
		pthread_mutexattr_t attr;
		int rc;

		rc = pthread_mutexattr_init(&attr);
		assert(rc == 0);

#ifdef PTHREAD_MUTEX_RECURSIVE
		(void)pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
#endif

		rc = pthread_mutex_init(&(ctx->mutex), &attr);
		assert(rc == 0);

		(void)pthread_mutexattr_destroy(&attr);
	}

	/* Initialize data */
	{
		ctx->tag = VALID_TAG_DEVICE;

		const qmckl_context_device context = (qmckl_context_device)ctx;
		assert(qmckl_context_check_device(context) !=
			   QMCKL_NULL_CONTEXT_DEVICE);

		qmckl_exit_code_device rc;

		ctx->numprec.precision = QMCKL_DEFAULT_PRECISION_DEVICE;
		ctx->numprec.range = QMCKL_DEFAULT_RANGE_DEVICE;

		rc = qmckl_init_point_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_electron_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_nucleus_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_ao_basis_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_mo_basis_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_determinant_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);

		rc = qmckl_init_jastrow_device(context);
		assert(rc == QMCKL_SUCCESS_DEVICE);
	}

	/* Allocate the host qmckl_memory_struct */
	{
		const size_t size = 128L;
		qmckl_memory_info_struct_device *new_array =
			calloc(size, sizeof(qmckl_memory_info_struct_device));
		if (new_array == NULL) {
			free(ctx);
			return QMCKL_NULL_CONTEXT_DEVICE;
		}
		memset(&(new_array[0]), 0,
			   size * sizeof(qmckl_memory_info_struct_device));

		ctx->memory.element = new_array;
		ctx->memory.array_size = size;
		ctx->memory.n_allocated = (size_t)0;
	}

	/* Allocate the device qmckl_memory_struct */
	const size_t size = 128L;
	qmckl_memory_info_struct_device *new_array =
		calloc(size, sizeof(qmckl_memory_info_struct_device));
	if (new_array == NULL) {
		free(ctx);
		return QMCKL_NULL_CONTEXT_DEVICE;
	}
	memset(&(new_array[0]), 0, size * sizeof(qmckl_memory_info_struct_device));

	ctx->memory_device.element = new_array;
	ctx->memory_device.array_size = size;
	ctx->memory_device.n_allocated = (size_t)0;

	/* Set the device_id */
	ctx->device_id = device_id;

	return (qmckl_context_device)ctx;
}

//**********
// CONTEXT DESTROY
//**********

qmckl_exit_code_device
qmckl_context_destroy_device(const qmckl_context_device context) {

	const qmckl_context_device checked_context =
		qmckl_context_check_device(context);
	if (checked_context == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_INVALID_CONTEXT_DEVICE;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx !=
		   NULL); /* Shouldn't be possible because the context is valid */

	int device_id = qmckl_get_device_id(context);

	/* Memory: Remove all allocated data */
	for (size_t pos = (size_t)0; pos < ctx->memory.array_size; ++pos) {
		if (ctx->memory.element[pos].pointer != NULL) {
			qmckl_free_device(context, ctx->memory.element[pos].pointer);
		}
	}
	assert(ctx->memory.n_allocated == (size_t)0);
	free(ctx->memory.element);
	ctx->memory.element = NULL;
	ctx->memory.array_size = (size_t)0;

	/* Device memory: Remove all allocated data */
	for (size_t pos = (size_t)0; pos < ctx->memory_device.array_size; ++pos) {
		if (ctx->memory_device.element[pos].pointer != NULL) {
			qmckl_free_device(context, ctx->memory_device.element[pos].pointer);
		}
	}
	assert(ctx->memory_device.n_allocated == (size_t)0);
	free(ctx->memory_device.element);
	ctx->memory_device.element = NULL;
	ctx->memory_device.array_size = (size_t)0;

	ctx->tag = INVALID_TAG_DEVICE;

	const int rc_destroy = pthread_mutex_destroy(&(ctx->mutex));
	if (rc_destroy != 0) {
		/* DEBUG */
		fprintf(stderr, "qmckl_context_destroy_device: %s (count = %d)\n",
				strerror(rc_destroy), ctx->lock_count);
		abort();
	}

	free(ctx);

	return QMCKL_SUCCESS_DEVICE;
}
