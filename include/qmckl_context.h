#pragma once

// This file contains prototypes for device context related functions,
// as well as the definition of qmckl_context_device and
// qmckl_context_device_struct, the device specific datatypes for context

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_memory.h"


qmckl_exit_code_device
qmckl_context_touch_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);
qmckl_exit_code_device
qmckl_context_destroy_device(const qmckl_context_device context);

static inline size_t qmckl_get_device_id(qmckl_context_device context) {
	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	return ctx->device_id;
}

qmckl_context_device
qmckl_context_check_device(const qmckl_context_device context);

void qmckl_lock_device(qmckl_context_device context);
void qmckl_unlock_device(qmckl_context_device context);
