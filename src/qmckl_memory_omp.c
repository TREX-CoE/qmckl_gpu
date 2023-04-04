#include "../include/qmckl_memory.h"
#include <assert.h>

// This file contains functions prototypes for context memory management
// functions (on device only, we expect most if not all of the context
// memory to be allocated on device in most cases)
// (OpenMP implementations)

//**********
// ALLOCS / FREES
//**********

void *qmckl_malloc_device(qmckl_context_device context, size_t size) {

	assert(qmckl_context_check_device(context != QMCKL_NULL_CONTEXT_DEVICE));

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	int device_id = qmckl_get_device_id(context);

	/* Allocate memory and zero it */
	void *pointer = omp_target_alloc(size, device_id);
	if (pointer == NULL) {
		return NULL;
	}

	// TODO
	// Memset to 0 of size info.size
	// memset(pointer, 0, info.size);

	qmckl_lock_device(context);
	{

		/* If qmckl_memory_struct is full, reallocate a larger one */
		if (ctx->memory_device.n_allocated == ctx->memory_device.array_size) {
			const size_t old_size = ctx->memory_device.array_size;
			qmckl_memory_info_struct_device *new_array =
				realloc(ctx->memory_device.element,
						2L * old_size * sizeof(qmckl_memory_info_struct_device));
			if (new_array == NULL) {
				qmckl_unlock_device(context);
				free(pointer);
				return NULL;
			}

			memset(&(new_array[old_size]), 0,
				   old_size * sizeof(qmckl_memory_info_struct_device));
			ctx->memory_device.element = new_array;
			ctx->memory_device.array_size = 2L * old_size;
		}

		/* Find first NULL entry */
		size_t pos = (size_t)0;
		while (pos < ctx->memory_device.array_size &&
			   ctx->memory_device.element[pos].size > (size_t)0) {
			pos += (size_t)1;
		}
		assert(ctx->memory_device.element[pos].size == (size_t)0);

		/* Copy info at the new location */
		ctx->memory_device.element[pos].size = size;
		ctx->memory_device.element[pos].pointer = pointer;
		ctx->memory_device.n_allocated += (size_t)1;
	}
	qmckl_unlock_device(context);

	return pointer;
}

qmckl_exit_code_device qmckl_free_device(qmckl_context_device context,
								  void *const ptr) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_free_device", NULL);
	}

	if (ptr == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
							  "qmckl_free_device", "NULL pointer");
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	int device_id = qmckl_get_device_id(context);

	qmckl_lock_device(context);
	{
		/* Find pointer in array of saved pointers */
		size_t pos = (size_t)0;
		while (pos < ctx->memory_device.array_size &&
			   ctx->memory_device.element[pos].pointer != ptr) {
			pos += (size_t)1;
		}

		if (pos >= ctx->memory_device.array_size) {
			/* Not found */
			qmckl_unlock_device(context);
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
								  "qmckl_free_device",
								  "Pointer not found in context");
		}

		omp_target_free(ptr, device_id);

		memset(&(ctx->memory_device.element[pos]), 0, sizeof(qmckl_memory_info_struct_device));
		ctx->memory_device.n_allocated -= (size_t)1;
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// MEMCPYS
//**********

qmckl_exit_code_device qmckl_memcpy_H2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_memcpy_H2D", NULL);
	}

	if (dest == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
							  "qmckl_memcpu_H2D", "NULL dest pointer");
	}

	if (src == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_memcpu_H2D", "NULL src pointer");
	}

	int device_id = qmckl_get_device_id(context);

	qmckl_lock_device(context);
	{
		int ret = omp_target_memcpy(dest, src, size, 0, 0, device_id,
									omp_get_initial_device());
		if (ret) {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
								  "qmckl_memcpy_H2D",
								  "Call to omp_target_memcpy failed");
		}
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_memcpy_D2H(qmckl_context_device context, void *const dest,
								 void *const src, size_t size) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_memcpy_D2H", NULL);
	}

	if (dest == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
							  "qmckl_memcpy_D2H", "NULL dest pointer");
	}

	if (src == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_memcpy_D2H", "NULL src pointer");
	}

	int device_id = qmckl_get_device_id(context);

	qmckl_lock_device(context);
	{
		int ret = omp_target_memcpy(dest, src, size, 0, 0,
									omp_get_initial_device(), device_id);
		if (ret) {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
								  "qmckl_memcpy_D2H",
								  "Call to omp_target_memcpy failed");
		}
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_memcpy_D2D(qmckl_context_device context, void *dest,
								 void *src, size_t size) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
							  "qmckl_memcpy_D2D", NULL);
	}

	if (dest == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
							  "qmckl_memcpy_D2D", "NULL dest pointer");
	}

	if (src == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
							  "qmckl_memcpy_D2D", "NULL src pointer");
	}

	int device_id = qmckl_get_device_id(context);

	qmckl_lock_device(context);
	{
		int ret =
			omp_target_memcpy(dest, src, size, 0, 0, device_id, device_id);
		if (ret) {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
								  "qmckl_memcpy_D2D",
								  "Call to omp_target_memcpy failed");
		}
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}
