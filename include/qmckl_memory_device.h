#pragma once

#include <qmckl.h>

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_memory_private_type.h"

#include "qmckl_context_device.h"
#include "qmckl_device_types.h"

void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr);

void *qmckl_malloc_omp_device(qmckl_context_device context,
							  const qmckl_memory_info_struct info,
							  int device_id);

qmckl_exit_code qmckl_free_omp_device(qmckl_context_device context,
									  void *const ptr, int device_id);
