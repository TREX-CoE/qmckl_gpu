#include "../include/qmckl_context_device.h"

// This file contains qmckl_context_device related functions


//**********
// CONTEXT DESTROY
//**********

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context) {

	const qmckl_context_device checked_context =
		qmckl_context_check((qmckl_context)context);
	if (checked_context == QMCKL_NULL_CONTEXT)
		return QMCKL_INVALID_CONTEXT;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx !=
		   NULL); /* Shouldn't be possible because the context is valid */

	qmckl_context_device_struct *const ds =
		(qmckl_context_device_struct *)ctx->qmckl_extra;
	int device_id = qmckl_get_device_id(context);

	qmckl_lock((qmckl_context)context);
	{
		/* Host memory: Remove all allocated data */
		for (size_t pos = (size_t)0; pos < ctx->memory.array_size; ++pos) {
			if (ctx->memory.element[pos].pointer != NULL) {
				omp_target_free(ctx->memory.element[pos].pointer, device_id);
				memset(&(ctx->memory.element[pos]), 0,
					   sizeof(qmckl_memory_info_struct));
				ctx->memory.n_allocated -= 1;
			}
		}
		assert(ctx->memory.n_allocated == (size_t)0);
		free(ctx->memory.element);
		ctx->memory.element = NULL;
		ctx->memory.array_size = (size_t)0;

		/* Device memory: Remove all allocated data */
		for (size_t pos = (size_t)0; pos < ds->memory.array_size; ++pos) {
			if (ds->memory.element[pos].pointer != NULL) {
				omp_target_free(ds->memory.element[pos].pointer, device_id);
				memset(&(ds->memory.element[pos]), 0,
					   sizeof(qmckl_memory_info_struct));
				ds->memory.n_allocated -= 1;
			}
		}
		assert(ds->memory.n_allocated == (size_t)0);
		free(ds->memory.element);
		ds->memory.element = NULL;
		ds->memory.array_size = (size_t)0;

		/* Free the qmckl_context_device_structured allocated in qmckl_extra */
		free(ctx->qmckl_extra);
	}
	qmckl_unlock((qmckl_context)context);

	ctx->tag = INVALID_TAG;

	const int rc_destroy = pthread_mutex_destroy(&(ctx->mutex));
	if (rc_destroy != 0) {
		/* DEBUG */
		fprintf(stderr, "qmckl_context_destroy: %s (count = %d)\n",
				strerror(rc_destroy), ctx->lock_count);
		abort();
	}

	free(ctx);

	return QMCKL_SUCCESS;
}
