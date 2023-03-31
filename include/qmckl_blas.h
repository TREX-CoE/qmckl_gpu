#pragma once

// This file contains functions prototypes for BLAS related functions
// (mostly manipulation of the device, matrix and tensor types)

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"

qmckl_vector_device qmckl_vector_alloc_device(qmckl_context_device context,
											  const int64_t size);

qmckl_exit_code_device
qmckl_vector_device_free_device(qmckl_context_device context,
								qmckl_vector_device *vector);

qmckl_exit_code_device qmckl_vector_device_of_double_device(
	const qmckl_context_device context, const double *target,
	const int64_t size_max, qmckl_vector_device *vector_out);

qmckl_matrix_device qmckl_matrix_alloc_device(qmckl_context_device context,
											  const int64_t size1,
											  const int64_t size2);

qmckl_exit_code_device
qmckl_matrix_device_free_device(qmckl_context_device context,
								qmckl_matrix_device *matrix);

qmckl_matrix_device qmckl_matrix_set_device(qmckl_matrix_device matrix,
											double value);

qmckl_exit_code_device qmckl_matrix_device_of_double_device(
	const qmckl_context_device context, const double *target,
	const int64_t size_max, qmckl_matrix_device *matrix_out);

qmckl_exit_code_device qmckl_transpose_device(qmckl_context_device context,
											  const qmckl_matrix_device A,
											  qmckl_matrix_device At);

qmckl_tensor_device qmckl_tensor_alloc_device(qmckl_context_device context,
											  const int64_t order,
											  const int64_t *size);

qmckl_exit_code_device
qmckl_tensor_device_free_device(qmckl_context_device context,
								qmckl_tensor_device *tensor);

qmckl_tensor_device qmckl_tensor_set_device(qmckl_tensor_device tensor,
											double value);
