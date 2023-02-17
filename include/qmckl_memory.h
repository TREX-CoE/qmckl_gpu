#pragma once

// This file contains functions prototypes for memory management functions

#include <qmckl.h>

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_memory_private_type.h"

#include "qmckl_context.h"

/* Allocs & frees */

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_device(qmckl_context_device context, size_t size);

qmckl_exit_code qmckl_free_device(qmckl_context_device context,
								  void *const ptr);

/* Memcpys */

qmckl_exit_code qmckl_memcpy_H2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code qmckl_memcpy_D2H(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
qmckl_exit_code qmckl_memcpy_D2D(qmckl_context_device context, void *const dest,
								 void *const src, size_t size);
