#include "../include/qmckl_blas.h"




#ifdef HAVE_CUBLAS
#include <cublas.h>
#elif HAVE_ROCBLAS 
#include "rocblas.h" 
#endif

// This file provides OpenMP implementations of BLAS functions (mostly
// initialization and manipulation of vector, matrix, ... types). All functions
// accept device pointers.
// (OpenMP/OpenACC independent functions only)

//**********
// VECTOR
//**********

qmckl_vector_device qmckl_vector_alloc_device(qmckl_context_device context,
											  const int64_t size) {
	/* Should always be true by contruction */
	assert(size > (int64_t)0);

	qmckl_vector_device result;
	result.size = size;

	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;
	result.data = qmckl_malloc_device(context, size * sizeof(double));

	if (result.data == NULL) {
		result.size = (int64_t)0;
	}

	return result;
}

qmckl_exit_code_device qmckl_vector_free_device(qmckl_context_device context,
												qmckl_vector_device *vector) {
	/* Always true */
	assert(vector->data != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_free_device(context, vector->data);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	vector->size = (int64_t)0;
	vector->data = NULL;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_vector_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_vector_device *vector_out) {

	int device_id = qmckl_get_device_id(context);

	qmckl_vector_device vector = *vector_out;
	/* Always true by construction */
	assert(((qmckl_context_device)context) != QMCKL_NULL_CONTEXT_DEVICE);

	if (vector.size == 0) {
		// This error is thrown
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_vector_device_of_double",
									 "Vector not allocated");
	}

	if (vector.size != size_max) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_vector_device_of_double",
									 "Wrong vector size");
	}

	omp_target_memcpy(vector.data, target, vector.size * sizeof(double), 0, 0,
					  device_id, device_id);

	*vector_out = vector;
	return QMCKL_SUCCESS_DEVICE;
}

//**********
// MATRIX
//**********

qmckl_matrix_device qmckl_matrix_alloc_device(qmckl_context_device context,
											  const int64_t size1,
											  const int64_t size2) {
	/* Should always be true by contruction */
	assert(size1 * size2 > (int64_t)0);

	qmckl_matrix_device result;

	result.size[0] = size1;
	result.size[1] = size2;

	result.data =
		(double *)qmckl_malloc_device(context, size1 * size2 * sizeof(double));

	if (result.data == NULL) {
		result.size[0] = (int64_t)0;
		result.size[1] = (int64_t)0;
	}
	double *data = result.data;

	return result;
}

qmckl_exit_code_device qmckl_matrix_free_device(qmckl_context_device context,
												qmckl_matrix_device *matrix) {
	/* Always true */
	assert(matrix->data != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_free_device(context, matrix->data);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}
	matrix->data = NULL;
	matrix->size[0] = (int64_t)0;
	matrix->size[1] = (int64_t)0;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_matrix_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_matrix_device *matrix_out) {

	// (assuming the matrix is already allocated)

	int device_id = qmckl_get_device_id(context);

	qmckl_matrix_device matrix = *matrix_out;
	/* Always true by construction */
	assert(((qmckl_context_device)context) != QMCKL_NULL_CONTEXT_DEVICE);

	if (matrix.size[0] * matrix.size[1] == 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_matrix_device_of_double_device",
									 "Matrix not allocated");
	}

	if (matrix.size[0] * matrix.size[1] > size_max) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_matrix_device_of_double_device",
									 "Wrong vector size");
	}

	omp_target_memcpy(matrix.data, target, size_max * sizeof(double), 0, 0,
					  device_id, device_id);

	*matrix_out = matrix;
	return QMCKL_SUCCESS_DEVICE;
}

//**********
// TENSOR
//**********

qmckl_tensor_device qmckl_tensor_alloc_device(qmckl_context_device context,
											  const int64_t order,
											  const int64_t *size) {
	/* Should always be true by construction */
	assert(order > 0);
	assert(order <= QMCKL_TENSOR_ORDER_MAX_DEVICE);
	assert(size != NULL);

	qmckl_tensor_device result;
	result.order = order;

	int64_t prod_size = (int64_t)1;
	for (int64_t i = 0; i < order; ++i) {
		assert(size[i] > (int64_t)0);
		result.size[i] = size[i];
		prod_size *= size[i];
	}

	result.data =
		(double *)qmckl_malloc_device(context, prod_size * sizeof(double));

	if (result.data == NULL) {
		memset(&result, 0, sizeof(qmckl_tensor_device));
	}

	return result;
}

qmckl_exit_code_device qmckl_tensor_free_device(qmckl_context_device context,
												qmckl_tensor_device *tensor) {
	/* Always true */
	assert(tensor->data != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_free_device(context, tensor->data);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	// TODO Memset to 0
	// memset(tensor, 0, sizeof(qmckl_tensor_device));

	return QMCKL_SUCCESS_DEVICE;
}





qmckl_exit_code_device qmckl_gpu_dgemm(qmckl_context_device context, char transA, char transB, int64_t m, int64_t n, int64_t k, double alpha, double* A, int64_t lda,\
 double* B, int64_t ldb, double beta, double* C, int64_t ldc)
{
    printf("coucou\n");

#pragma omp target is_device_ptr(A, B, C)
{
    #ifdef HAVE_CUBLAS
    cublasStatus ret = CUBLAS_STATUS_SUCCESS;
    cublasHandle_t handle;
    ret = cublasDgemm(handle, transA, transB, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
    if( ret != CUBLAS_STATUS_SUCCES)
    {
        retrun  QMCKL_FAILURE_DEVICE;
    }

    #elif HAVE_ROCBLAS
    rocblas_status ret = rocblas_status_success;
    rocblas_handle = handle;
    rocblas_create_handle(&handle);
    ret = rocblas_dgemm(handle, transA, transB, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
    if(ret != rocblas_status_succes)
    {
        return QMCKL_FAILURE_DEVICE;
    }
    #else
    #pragma omp teams distribute parallel for
    for(int i = 0; i < lda; i++)
    {
        for(int j = 0; j < ldb; j++)
        {
            C[i * lda + j] = 0.0;
            for(int k = 0; k < ldc; k++)
            {
                C[i * lda + j] += A[i * lda + k] * B[j * ldb + k];
            }
        }
    }
    #endif
}
    return QMCKL_SUCCESS_DEVICE;    
}