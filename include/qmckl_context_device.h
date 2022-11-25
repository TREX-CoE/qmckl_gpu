#pragma once

// This file contains prototypes for device context related functions,
// ass well as the definition of qmckl_context_device_struct, which
// contains additional and device-related data.

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

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_blas_device.h"
#include "qmckl_device_types.h"
#include "qmckl_memory_device.h"

typedef struct {
	size_t device_id;
	qmckl_memory_struct memory;
} qmckl_context_device_struct;

qmckl_exit_code
qmckl_context_destroy_omp_device(const qmckl_context_device context);

qmckl_exit_code qmckl_context_touch_device(const qmckl_context_device context);
qmckl_exit_code
qmckl_context_touch_omp_device(const qmckl_context_device context);
qmckl_exit_code
qmckl_context_touch_acc_device(const qmckl_context_device context);

qmckl_context_device qmckl_context_create_device();
qmckl_context_device qmckl_context_create_omp_device();
qmckl_context_device qmckl_context_create_acc_device();

static inline size_t qmckl_get_device_id(qmckl_context_device context) {
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	qmckl_context_device_struct *const ds =
		(qmckl_context_device_struct *)ctx->qmckl_extra;
	return ds->device_id;
}
