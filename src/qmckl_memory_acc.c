#include "../include/qmckl_memory.h"
#include <assert.h>
#include <openacc.h>

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
	void *pointer = acc_malloc(size);
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

		acc_free(ptr);

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
	{ acc_memcpy_to_device(dest, src, size); }
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
	{ acc_memcpy_from_device(dest, src, size); }
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
// NOT working from device to device
// acc_memcpy_to_device(dest, src, size);
//
// NOT supporting device-pointers only
// int dest_dev_id = 0;
// int src_dev_id = 0;
// acc_memcpy_d2d(dest, src, size, dest_dev_id, src_dev_id);
//
// Workaround (cast as int if multiple of 4 bytes?)
#pragma acc parallel loop gang vector deviceptr(src, dest), copyin(size)
		for (int i = 0; i < size; ++i)
			((char *)dest)[i] = ((char *)src)[i];
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}
