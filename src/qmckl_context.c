// This file contains qmckl_context_device related functions

#include "../include/qmckl_context.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code qmckl_context_touch_device(const qmckl_context_device context) {
	return qmckl_context_touch((qmckl_context)context);
}

//**********
// CONTEXT CREATE
//**********

qmckl_context_device qmckl_context_create_device(int device_id) {
	qmckl_context_device context = (qmckl_context_device)qmckl_context_create();
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	/* Allocate the qmckl_context_device_struct */
	ctx->qmckl_extra = malloc(sizeof(qmckl_context_device_struct));

	qmckl_context_device_struct *const ds =
		(qmckl_context_device_struct *)ctx->qmckl_extra;

	/* Allocate the device qmckl_memory_struct */
	const size_t size = 128L;
	qmckl_memory_info_struct *new_array =
		calloc(size, sizeof(qmckl_memory_info_struct));
	if (new_array == NULL) {
		free(ctx);
		return QMCKL_NULL_CONTEXT;
	}
	memset(&(new_array[0]), 0, size * sizeof(qmckl_memory_info_struct));

	ds->memory.element = new_array;
	ds->memory.array_size = size;
	ds->memory.n_allocated = (size_t)0;

	/* Set the device_id */
	ds->device_id = device_id;

	return context;
}

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
		free(ctx->qmckl_extra);
		/* Memory: Remove all allocated data */
		for (size_t pos = (size_t)0; pos < ctx->memory.array_size; ++pos) {
			if (ctx->memory.element[pos].pointer != NULL) {
				qmckl_free(context, ctx->memory.element[pos].pointer);
			}
		}
		assert(ctx->memory.n_allocated == (size_t)0);
		free(ctx->memory.element);
		ctx->memory.element = NULL;
		ctx->memory.array_size = (size_t)0;

		/* Device memory: Remove all allocated data */
		for (size_t pos = (size_t)0; pos < ds->memory.array_size; ++pos) {
			if (ds->memory.element[pos].pointer != NULL) {
				qmckl_free_device(context, ds->memory.element[pos].pointer);
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
