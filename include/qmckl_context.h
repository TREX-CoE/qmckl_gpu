#pragma once

// This file contains prototypes for device context related functions,
// as well as the definition of qmckl_context_device and
// qmckl_context_device_struct, the device specific datatypes for context

#include <qmckl.h>

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "qmckl_basic_functions.h"

typedef int64_t qmckl_context_device;

#include "qmckl_memory.h"

typedef struct {
	size_t device_id;
	qmckl_memory_struct memory;
} qmckl_context_device_struct;

qmckl_exit_code
qmckl_context_destroy_omp_device(const qmckl_context_device context);

qmckl_exit_code qmckl_context_touch_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);

static inline size_t qmckl_get_device_id(qmckl_context_device context) {
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	qmckl_context_device_struct *const ds =
		(qmckl_context_device_struct *)ctx->qmckl_extra;
	return ds->device_id;
}
