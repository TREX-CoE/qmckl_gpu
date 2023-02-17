#include "../include/qmckl_blas.h"

// This file provides OpenMP implementations of BLAS functions (mostly
// initialization and manipulation of vector, matrix, ... types). All functions
// accept device pointers.
// (OpenMP/OpenACC independent functions only)

//**********
// VECTOR
//**********

qmckl_vector qmckl_vector_alloc_device(qmckl_context_device context,
									   const int64_t size) {
	/* Should always be true by contruction */
	assert(size > (int64_t)0);

	qmckl_vector result;
	result.size = size;

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	result.data = qmckl_malloc_device(context, size * sizeof(double));

	if (result.data == NULL) {
		result.size = (int64_t)0;
	}

	return result;
}

qmckl_exit_code qmckl_vector_free_device(qmckl_context_device context,
										 qmckl_vector *vector) {
	/* Always true */
	assert(vector->data != NULL);

	qmckl_exit_code rc;

	rc = qmckl_free_device(context, vector->data);
	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	vector->size = (int64_t)0;
	vector->data = NULL;
	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_vector_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_vector *vector_out) {

	int device_id = qmckl_get_device_id(context);

	qmckl_vector vector = *vector_out;
	/* Always true by construction */
	assert(qmckl_context_check((qmckl_context)context) != QMCKL_NULL_CONTEXT);

	if (vector.size == 0) {
		// This error is thrown
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_vector_of_double", "Vector not allocated");
	}

	if (vector.size != size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_vector_of_double", "Wrong vector size");
	}

	omp_target_memcpy(vector.data, target, vector.size * sizeof(double), 0, 0,
					  device_id, device_id);

	*vector_out = vector;
	return QMCKL_SUCCESS;
}

//**********
// MATRIX
//**********

qmckl_matrix qmckl_matrix_alloc_device(qmckl_context_device context,
									   const int64_t size1,
									   const int64_t size2) {
	/* Should always be true by contruction */
	assert(size1 * size2 > (int64_t)0);

	qmckl_matrix result;

	result.size[0] = size1;
	result.size[1] = size2;

	result.data = (double *)qmckl_malloc_device(context, size1 * size2 * sizeof(double));

	if (result.data == NULL) {
		result.size[0] = (int64_t)0;
		result.size[1] = (int64_t)0;
	}

	return result;
}

qmckl_exit_code qmckl_matrix_free_device(qmckl_context_device context,
										 qmckl_matrix *matrix) {
	/* Always true */
	assert(matrix->data != NULL);

	qmckl_exit_code rc;

	rc = qmckl_free_device(context, matrix->data);
	if (rc != QMCKL_SUCCESS) {
		return rc;
	}
	matrix->data = NULL;
	matrix->size[0] = (int64_t)0;
	matrix->size[1] = (int64_t)0;

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_matrix_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_matrix *matrix_out) {

	// (assuming the matrix is already allocated)

	int device_id = qmckl_get_device_id(context);

	qmckl_matrix matrix = *matrix_out;
	/* Always true by construction */
	assert(qmckl_context_check((qmckl_context)context) != QMCKL_NULL_CONTEXT);

	if (matrix.size[0] * matrix.size[1] == 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_matrix_of_double_device",
							  "Matrix not allocated");
	}

	if (matrix.size[0] * matrix.size[1] > size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_matrix_of_double_device",
							  "Wrong vector size");
	}

	omp_target_memcpy(matrix.data, target, size_max * sizeof(double), 0, 0,
					  device_id, device_id);

	*matrix_out = matrix;
	return QMCKL_SUCCESS;
}

//**********
// TENSOR
//**********

qmckl_tensor qmckl_tensor_alloc_device(qmckl_context context,
									   const int64_t order,
									   const int64_t *size) {
	/* Should always be true by construction */
	assert(order > 0);
	assert(order <= QMCKL_TENSOR_ORDER_MAX);
	assert(size != NULL);

	qmckl_tensor result;
	result.order = order;

	int64_t prod_size = (int64_t)1;
	for (int64_t i = 0; i < order; ++i) {
		assert(size[i] > (int64_t)0);
		result.size[i] = size[i];
		prod_size *= size[i];
	}

	result.data = (double *)qmckl_malloc_device(context, prod_size * sizeof(double));

	if (result.data == NULL) {
		memset(&result, 0, sizeof(qmckl_tensor));
	}

	return result;
}

qmckl_exit_code qmckl_tensor_free_device(qmckl_context_device context,
										 qmckl_tensor *tensor) {
	/* Always true */
	assert(tensor->data != NULL);

	qmckl_exit_code rc;

	rc = qmckl_free_device(context, tensor->data);
	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	// TODO Memset to 0
	// memset(tensor, 0, sizeof(qmckl_tensor));

	return QMCKL_SUCCESS;
}
