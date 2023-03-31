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

// TODO Replace by file content
#include "qmckl_basic_type.h"

typedef struct qmckl_context_struct_device {
	/* -- State of the library -- */

	/* Device id, only used w/ OpenMP */
	size_t device_id;

	/* Validity checking */
	uint64_t tag;

	/* Numerical precision */
	qmckl_numprec_struct_device numprec;

	/* Thread lock */
	int lock_count;
	pthread_mutex_t mutex;

	/* Error handling */
	qmckl_error_struct_device error;

	/* Memory allocation */
	qmckl_memory_struct_device memory;
	qmckl_memory_struct_device memory_device;

	/* Current date */
	uint64_t date;

	/* Points */
	qmckl_point_struct_device point;

	/* -- Molecular system -- */
	qmckl_nucleus_struct_device nucleus;
	qmckl_electron_struct_device electron;
	qmckl_ao_basis_struct_device ao_basis;
	qmckl_mo_basis_struct_device mo_basis;
	// TODO When available, add the Jastrow struct
	// qmckl_jastrow_struct_device jastrow;
	qmckl_determinant_struct_device det;
	qmckl_local_energy_struct_device local_energy;

	/* To be implemented:
	 */

	/* Pointer to implementation-specific data */

	void *qmckl_extra;

} qmckl_context_struct_device;

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
