// This file will contains general functions accepting a qmckl_contect_device as
// argument, hence why they need a device alternative

#include "../include/qmckl_context_device.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code qmckl_context_touch_device(const qmckl_context_device context) {
	return qmckl_context_touch((qmckl_context)context);
}

//**********
// CONTEXT CREATE
//**********
// NOTE The context destroy function frees some memory, so its implementation is
// OpenMP/OpenACC dependent

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
