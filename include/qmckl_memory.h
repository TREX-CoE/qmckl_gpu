#pragma once

// This file contains functions prototypes for memory management functions

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "qmckl_types.h"
#include "qmckl_context.h"

/* Allocs & frees */
void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct_device info);

void *qmckl_malloc_device(qmckl_context_device context, int64_t size);

qmckl_exit_code_device qmckl_free_host(qmckl_context_device context,
									   void *const ptr);

qmckl_exit_code_device qmckl_free_device(qmckl_context_device context,
										 void *const ptr);

/* Memcpys */

qmckl_exit_code_device qmckl_memcpy_H2D(qmckl_context_device context,
										void *const dest, void *const src,
										int64_t size);


qmckl_exit_code_device qmckl_memcpy_H2D_double(qmckl_context_device context,
										double *const dest, double *const src,
										int64_t numberElements);


qmckl_exit_code_device qmckl_memcpy_H2D_int64(qmckl_context_device context,
										int64_t *const dest, double *const src,
										int64_t numberElements);


qmckl_exit_code_device qmckl_memcpy_H2D_int32(qmckl_context_device context,
										int32_t *const dest, double *const src,
										int64_t numberElements);


qmckl_exit_code_device qmckl_memcpy_D2H(qmckl_context_device context,
										void *const dest, void *const src,
										int64_t size);





qmckl_exit_code_device qmckl_memcpy_D2H_double(qmckl_context_device context,
										double *const dest, double *const src,
										int64_t numberElements);


qmckl_exit_code_device qmckl_memcpy_D2H_int64(qmckl_context_device context,
										int64_t *const dest, double *const src,
										int64_t numberElements);


qmckl_exit_code_device qmckl_memcpy_D2H_int32(qmckl_context_device context,
										int32_t *const dest, double *const src,
										int64_t numberElements);






qmckl_exit_code_device qmckl_memcpy_D2D(qmckl_context_device context,
										void *const dest, void *const src,
										size_t size);
