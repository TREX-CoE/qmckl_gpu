// This file will contains general functions accepting a qmckl_contect_device as
// argument, hence why they need a device alternative

#include "../include/qmckl_context_device.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code
qmckl_context_touch_omp_device(const qmckl_context_device context) {
	return qmckl_context_touch((qmckl_context)context);
}

//**********
// CONTEXT CREATE
//**********

qmckl_context_device qmckl_context_create_omp_device() {
	return (qmckl_context_device)qmckl_context_create();
}

//**********
// CONTEXT DESTROY
//**********

qmckl_exit_code
qmckl_context_destroy_omp_device(const qmckl_context_device context,
				 int device_id) {

	const qmckl_context_device checked_context =
	    qmckl_context_check((qmckl_context)context);
	if (checked_context == QMCKL_NULL_CONTEXT)
		return QMCKL_INVALID_CONTEXT;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx !=
	       NULL); /* Shouldn't be possible because the context is valid */

	qmckl_lock((qmckl_context)context);
	{
		/* Memory: Remove all allocated data */
		for (size_t pos = (size_t)0; pos < ctx->memory.array_size;
		     ++pos) {
			if (ctx->memory.element[pos].pointer != NULL) {
				omp_target_free(
				    ctx->memory.element[pos].pointer,
				    device_id);
				memset(&(ctx->memory.element[pos]), 0,
				       sizeof(qmckl_memory_info_struct));
				ctx->memory.n_allocated -= 1;
			}
		}
		assert(ctx->memory.n_allocated == (size_t)0);
		free(ctx->memory.element);
		ctx->memory.element = NULL;
		ctx->memory.array_size = (size_t)0;
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
