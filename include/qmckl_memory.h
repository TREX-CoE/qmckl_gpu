#pragma once

// This file contains functions prototypes for memory management functions

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "qmckl_basic_functions.h"

typedef int64_t qmckl_context_device;

typedef struct qmckl_memory_info_struct {
	size_t size;
	void* pointer;
} qmckl_memory_info_struct_device;

typedef struct qmckl_memory_struct_device {
	size_t n_allocated;
	size_t array_size;
	qmckl_memory_info_struct_device *element;
} qmckl_memory_struct_device;

static const qmckl_memory_info_struct_device qmckl_memory_info_struct_zero_device = {
	.size = (size_t)0,
	.pointer = NULL
};


/* Allocs & frees */

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct_device info);

qmckl_exit_code_device qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct_device info);

qmckl_exit_code_device qmckl_free_device(qmckl_context_device context,
								  void *const ptr);

/* Memcpys */

qmckl_exit_code_device qmckl_memcpy_H2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code_device qmckl_memcpy_D2H(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code_device qmckl_memcpy_D2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
