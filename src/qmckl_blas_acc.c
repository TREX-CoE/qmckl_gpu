#include "../include/qmckl_blas.h"
#include <acc.h>

// This file provides OpenMP implementations of BLAS functions (mostly
// initialization and manipulation of vector, matrix, ... types). All functions
// accept device pointers.
// (functions requiring OpenACC pragmas only)

//**********
// MATRIX
//**********

qmckl_matrix qmckl_matrix_set_device(qmckl_matrix matrix, double value) {
	// Recompute array size
	int prod_size = matrix.size[0] * matrix.size[1];

	double *data = matrix.data;
#pragma acc data deviceptr(data)
	{
		for (int i = 0; i < prod_size; i++) {
			data[i] = value;
		}
	}
	return matrix;
}

qmckl_exit_code qmckl_transpose_device(qmckl_context_device context,
									   const qmckl_matrix A, qmckl_matrix At) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (A.size[0] < 1) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_transpose_device", "Invalid size for A");
	}

	if (At.data == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_transpose_device",
							  "Output matrix not allocated");
	}

	if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_transpose_device", "Invalid size for At");
	}

	double *A_data = A.data;
	int A_s0 = A.size[0];

	double *At_data = At.data;
	int At_s0 = At.size[0];
	int At_s1 = At.size[1];

#pragma acc data deviceptr(A_data, At_data)
	{
#pragma acc parallel loop collapse(2)
		for (int64_t j = 0; j < At_s1; ++j)
			for (int64_t i = 0; i < At_s0; ++i)
				At_data[i + j * At_s0] = A_data[j + i * A_s0];
	}

	return QMCKL_SUCCESS;
}

//**********
// TENSOR
//**********

qmckl_tensor qmckl_tensor_set_device(qmckl_tensor tensor, double value) {
	// Recompute array size
	int prod_size = 1;

	for (int i = 0; i < tensor.order; i++) {
		prod_size *= tensor.size[i];
	}

	double *data = tensor.data;
#pragma acc data deviceptr(data)
	{
#pragma acc parallel loop
		for (int i = 0; i < prod_size; i++) {
			data[i] = value;
		}
	}
	return tensor;
}
