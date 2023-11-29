#include "../include/qmckl_mo.h"
#include "include/qmckl_ao.h"
#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#endif
#ifdef HAVE_CUSPARSE
#include <cuda_runtime_api.h>
#include <cusparse_v2.h>
#endif

/* Provided check  */

bool qmckl_mo_basis_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->mo_basis.provided;
}

//**********
// COMPUTE
//**********

#ifdef HAVE_CUBLAS
qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_sgemm_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_vgl,
	double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	cublasHandle_t handle;

	float const alpha = 1.;
	float const beta = 0.;

	cublasOperation_t transa = CUBLAS_OP_N;
	cublasOperation_t transb = CUBLAS_OP_N;

	// NB: cublas views arrays as column-major (ie FORTRAN-like) ordered
	int const m = mo_num;
	int const k = ao_num;
	int const n = point_num * 5;

	int const lda = m;
	int const ldb = k;
	int const ldc = m;

	float *A = qmckl_malloc_device(context, sizeof(float) * ao_num * mo_num);
	float *B =
		qmckl_malloc_device(context, sizeof(float) * 5 * ao_num * point_num);
	float *C =
		qmckl_malloc_device(context, sizeof(float) * 5 * mo_num * point_num);

#pragma omp target teams distribute parallel for simd is_device_ptr(           \
	A, coefficient_t) map(to                                                   \
						  : ao_num, mo_num)
#pragma acc parallel loop gang vector deviceptr(A, coefficient_t)              \
	copyin(ao_num, mo_num)
	for (int ii = 0; ii < ao_num * mo_num; ++ii) {
		A[ii] = (float)coefficient_t[ii];
	}
#pragma omp target teams distribute parallel for simd is_device_ptr(B, ao_vgl) \
	map(to                                                                     \
		: ao_num, point_num)
#pragma acc parallel loop gang vector deviceptr(B, ao_vgl)                     \
	copyin(ao_num, point_num)
	for (int ii = 0; ii < 5 * ao_num * point_num; ++ii) {
		B[ii] = (float)ao_vgl[ii];
	}

	cublasCreate(&handle);
	cublasSgemm_v2(handle, transa, transb, m, n, k, &alpha, A, lda, B, ldb,
				   &beta, C, ldc);
	cublasDestroy(handle);

#pragma omp target teams distribute parallel for simd is_device_ptr(C, mo_vgl) \
	map(to                                                                     \
		: mo_num, point_num)
#pragma acc parallel loop gang vector deviceptr(C, mo_vgl)                     \
	copyin(mo_num, point_num)
	for (int ii = 0; ii < 5 * mo_num * point_num; ++ii) {
		mo_vgl[ii] = (double)C[ii];
	}

	qmckl_free_device(context, A);
	qmckl_free_device(context, B);
	qmckl_free_device(context, C);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_dgemm_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_vgl,
	double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	cublasHandle_t handle;

	double const alpha = 1.;
	double const beta = 0.;

	cublasOperation_t transa = CUBLAS_OP_N;
	cublasOperation_t transb = CUBLAS_OP_N;

	// NB: cublas views arrays as column-major (ie FORTRAN-like) ordered
	int const m = mo_num;
	int const k = ao_num;
	int const n = point_num * 5;

	int const lda = m;
	int const ldb = k;
	int const ldc = m;

	cublasCreate(&handle);
	cublasDgemm_v2(handle, transa, transb, m, n, k, &alpha, coefficient_t, lda,
				   ao_vgl, ldb, &beta, mo_vgl, ldc);
	cublasDestroy(handle);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_dgemm_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_value,
	double *restrict const mo_value) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	cublasHandle_t handle;

	double const alpha = 1.;
	double const beta = 0.;

	cublasOperation_t transa = CUBLAS_OP_N;
	cublasOperation_t transb = CUBLAS_OP_N;

	// NB: cublas views arrays as column-major (ie FORTRAN-like) ordered
	int const m = mo_num;
	int const k = ao_num;
	int const n = point_num;

	int const lda = m;
	int const ldb = k;
	int const ldc = m;

	cublasCreate(&handle);
	cublasDgemm_v2(handle, transa, transb, m, n, k, &alpha, coefficient_t, lda,
				   ao_value, ldb, &beta, mo_value, ldc);
	cublasDestroy(handle);

	return QMCKL_SUCCESS_DEVICE;
}
qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_sgemm_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_vgl,
	double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	cublasHandle_t handle;

	float const alpha = 1.;
	float const beta = 0.;

	cublasOperation_t transa = CUBLAS_OP_N;
	cublasOperation_t transb = CUBLAS_OP_N;

	// NB: cublas views arrays as column-major (ie FORTRAN-like) ordered
	int const m = mo_num;
	int const k = ao_num;
	int const n = point_num;

	int const lda = m;
	int const ldb = k;
	int const ldc = m;

	float *A = qmckl_malloc_device(context, sizeof(float) * ao_num * mo_num);
	float *B = qmckl_malloc_device(context, sizeof(float) * ao_num * point_num);
	float *C = qmckl_malloc_device(context, sizeof(float) * mo_num * point_num);

#pragma omp target teams distribute parallel for simd is_device_ptr(           \
	A, coefficient_t) map(to                                                   \
						  : ao_num, mo_num)
#pragma acc parallel loop gang vector deviceptr(A, coefficient_t)              \
	copyin(ao_num, mo_num)
	for (int ii = 0; ii < ao_num * mo_num; ++ii) {
		A[ii] = (float)coefficient_t[ii];
	}
#pragma omp target teams distribute parallel for simd is_device_ptr(B, ao_vgl) \
	map(to                                                                     \
		: ao_num, point_num)
#pragma acc parallel loop gang vector deviceptr(B, ao_vgl)                     \
	copyin(ao_num, point_num)
	for (int ii = 0; ii < ao_num * point_num; ++ii) {
		B[ii] = (float)ao_vgl[ii];
	}

	cublasCreate(&handle);
	cublasSgemm_v2(handle, transa, transb, m, n, k, &alpha, A, lda, B, ldb,
				   &beta, C, ldc);
	cublasDestroy(handle);

#pragma omp target teams distribute parallel for simd is_device_ptr(C, mo_vgl) \
	map(to                                                                     \
		: mo_num, point_num)
#pragma acc parallel loop gang vector deviceptr(C, mo_vgl)                     \
	copyin(mo_num, point_num)
	for (int ii = 0; ii < mo_num * point_num; ++ii) {
		mo_vgl[ii] = (double)C[ii];
	}

	qmckl_free_device(context, A);
	qmckl_free_device(context, B);
	qmckl_free_device(context, C);

	return QMCKL_SUCCESS_DEVICE;
}
#endif

#ifdef HAVE_CUSPARSE
qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_cusparse_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_vgl,
	double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	// cusparse (dense) matrix descriptors
	//
	cusparseDnMatDescr_t matAd;
	int64_t A_nrows = point_num * 5;
	int64_t A_ncols = ao_num;
	int64_t lda = A_ncols;
	cusparseOrder_t A_storage_order = CUSPARSE_ORDER_ROW;

	cusparseDnMatDescr_t matBd;
	int64_t B_nrows = ao_num;
	int64_t B_ncols = mo_num;
	int64_t ldb = B_ncols;
	cusparseOrder_t B_storage_order = CUSPARSE_ORDER_ROW;

	cusparseDnMatDescr_t matCd;
	int64_t C_nrows = point_num * 5;
	int64_t C_ncols = mo_num;
	int64_t ldc = C_ncols;
	cusparseOrder_t C_storage_order = CUSPARSE_ORDER_ROW;

	cudaDataType cuda_datatype = CUDA_R_64F;

	cusparseHandle_t handle;
	cusparseCreate(&handle);

	cusparseCreateDnMat(&matAd, A_nrows, A_ncols, lda, ao_vgl, cuda_datatype,
						A_storage_order);
	cusparseCreateDnMat(&matBd, B_nrows, B_ncols, ldb, coefficient_t,
						cuda_datatype, B_storage_order);
	cusparseCreateDnMat(&matCd, C_nrows, C_ncols, ldc, mo_vgl, cuda_datatype,
						C_storage_order);

	// Convert sparse matrix A from dense to sparse(csr) format
	//
	cusparseSpMatDescr_t matAs;
	int64_t As_nnz = 0;
	int64_t *As_csr_offsets = NULL;
	int64_t *As_csr_columns = NULL;
	double *As_csr_values = NULL;

	cusparseIndexBase_t As_csrRowOffsetsType = CUSPARSE_INDEX_64I;
	cusparseIndexBase_t As_csrColIndType = CUSPARSE_INDEX_64I;
	cusparseIndexBase_t As_idxBase = CUSPARSE_INDEX_BASE_ZERO;

	// Why allocations needs to be done here?
	cudaMalloc((void **)&As_csr_offsets, (A_nrows + 1) * sizeof(int64_t));

	cusparseCreateCsr(&matAs, A_nrows, A_ncols, As_nnz, As_csr_offsets,
					  As_csr_columns, As_csr_values, As_csrRowOffsetsType,
					  As_csrColIndType, As_idxBase, cuda_datatype);

	size_t As_buffer_size;
	void *As_buffer = NULL;
	cusparseDenseToSparseAlg_t AlgDtoS = CUSPARSE_DENSETOSPARSE_ALG_DEFAULT;

	cusparseDenseToSparse_bufferSize(handle, matAd, matAs, AlgDtoS,
									 &As_buffer_size);
	cudaMalloc(&As_buffer, As_buffer_size);

	cusparseDenseToSparse_analysis(handle, matAd, matAs, AlgDtoS, As_buffer);
	cusparseSpMatGetSize(matAs, &A_nrows, &A_ncols, &As_nnz);

	cudaMalloc((void **)&As_csr_columns, As_nnz * sizeof(int64_t));
	cudaMalloc((void **)&As_csr_values, As_nnz * sizeof(double));

	cusparseCsrSetPointers(matAs, As_csr_offsets, As_csr_columns,
						   As_csr_values);
	cusparseDenseToSparse_convert(handle, matAd, matAs, AlgDtoS, As_buffer);

	// Sparse A * dense B = dense C
	//
	double const alpha = 1.;
	double const beta = 0.;
	cusparseOperation_t opA = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseOperation_t opB = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseSpMMAlg_t AlgSpMM = CUSPARSE_SPMM_ALG_DEFAULT;

	size_t SpMM_buffer_size;
	void *SpMM_buffer = NULL;

	cusparseSpMM_bufferSize(handle, opA, opB, &alpha, matAs, matBd, &beta,
							matCd, cuda_datatype, AlgSpMM, &SpMM_buffer_size);
	cudaMalloc(&SpMM_buffer, SpMM_buffer_size);

	cusparseSpMM(handle, opA, opB, &alpha, matAs, matBd, &beta, matCd,
				 cuda_datatype, AlgSpMM, SpMM_buffer);

	cudaFree(As_csr_offsets);
	cudaFree(As_buffer);
	cudaFree(As_csr_columns);
	cudaFree(As_csr_values);
	cudaFree(SpMM_buffer);

	cusparseDestroyDnMat(matAd);
	cusparseDestroyDnMat(matBd);
	cusparseDestroyDnMat(matCd);
	cusparseDestroySpMat(matAs);
	cusparseDestroy(handle);

	return QMCKL_SUCCESS_DEVICE;
}
qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_cusparse_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t mo_num, const int64_t point_num,
	const double *restrict coefficient_t, const double *restrict ao_value,
	double *restrict const mo_value) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	// cusparse (dense) matrix descriptors
	//
	cusparseDnMatDescr_t matAd;
	int64_t A_nrows = point_num;
	int64_t A_ncols = ao_num;
	int64_t lda = A_ncols;
	cusparseOrder_t A_storage_order = CUSPARSE_ORDER_ROW;

	cusparseDnMatDescr_t matBd;
	int64_t B_nrows = ao_num;
	int64_t B_ncols = mo_num;
	int64_t ldb = B_ncols;
	cusparseOrder_t B_storage_order = CUSPARSE_ORDER_ROW;

	cusparseDnMatDescr_t matCd;
	int64_t C_nrows = point_num;
	int64_t C_ncols = mo_num;
	int64_t ldc = C_ncols;
	cusparseOrder_t C_storage_order = CUSPARSE_ORDER_ROW;

	cudaDataType cuda_datatype = CUDA_R_64F;

	cusparseHandle_t handle;
	cusparseCreate(&handle);

	cusparseCreateDnMat(&matAd, A_nrows, A_ncols, lda, ao_value, cuda_datatype,
						A_storage_order);
	cusparseCreateDnMat(&matBd, B_nrows, B_ncols, ldb, coefficient_t,
						cuda_datatype, B_storage_order);
	cusparseCreateDnMat(&matCd, C_nrows, C_ncols, ldc, mo_value, cuda_datatype,
						C_storage_order);

	// Convert sparse matrix A from dense to sparse(csr) format
	//
	cusparseSpMatDescr_t matAs;
	int64_t As_nnz = 0;
	int64_t *As_csr_offsets = NULL;
	int64_t *As_csr_columns = NULL;
	double *As_csr_values = NULL;

	cusparseIndexBase_t As_csrRowOffsetsType = CUSPARSE_INDEX_64I;
	cusparseIndexBase_t As_csrColIndType = CUSPARSE_INDEX_64I;
	cusparseIndexBase_t As_idxBase = CUSPARSE_INDEX_BASE_ZERO;

	// Why allocations needs to be done here?
	cudaMalloc((void **)&As_csr_offsets, (A_nrows + 1) * sizeof(int64_t));

	cusparseCreateCsr(&matAs, A_nrows, A_ncols, As_nnz, As_csr_offsets,
					  As_csr_columns, As_csr_values, As_csrRowOffsetsType,
					  As_csrColIndType, As_idxBase, cuda_datatype);

	size_t As_buffer_size;
	void *As_buffer = NULL;
	cusparseDenseToSparseAlg_t AlgDtoS = CUSPARSE_DENSETOSPARSE_ALG_DEFAULT;

	cusparseDenseToSparse_bufferSize(handle, matAd, matAs, AlgDtoS,
									 &As_buffer_size);
	cudaMalloc(&As_buffer, As_buffer_size);

	cusparseDenseToSparse_analysis(handle, matAd, matAs, AlgDtoS, As_buffer);
	cusparseSpMatGetSize(matAs, &A_nrows, &A_ncols, &As_nnz);

	cudaMalloc((void **)&As_csr_columns, As_nnz * sizeof(int64_t));
	cudaMalloc((void **)&As_csr_values, As_nnz * sizeof(double));

	cusparseCsrSetPointers(matAs, As_csr_offsets, As_csr_columns,
						   As_csr_values);
	cusparseDenseToSparse_convert(handle, matAd, matAs, AlgDtoS, As_buffer);

	// Sparse A * dense B = dense C
	//
	double const alpha = 1.;
	double const beta = 0.;
	cusparseOperation_t opA = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseOperation_t opB = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseSpMMAlg_t AlgSpMM = CUSPARSE_SPMM_ALG_DEFAULT;

	size_t SpMM_buffer_size;
	void *SpMM_buffer = NULL;

	cusparseSpMM_bufferSize(handle, opA, opB, &alpha, matAs, matBd, &beta,
							matCd, cuda_datatype, AlgSpMM, &SpMM_buffer_size);
	cudaMalloc(&SpMM_buffer, SpMM_buffer_size);

	cusparseSpMM(handle, opA, opB, &alpha, matAs, matBd, &beta, matCd,
				 cuda_datatype, AlgSpMM, SpMM_buffer);

	cudaFree(As_csr_offsets);
	cudaFree(As_buffer);
	cudaFree(As_csr_columns);
	cudaFree(As_csr_values);
	cudaFree(SpMM_buffer);

	cusparseDestroyDnMat(matAd);
	cusparseDestroyDnMat(matBd);
	cusparseDestroyDnMat(matCd);
	cusparseDestroySpMat(matAs);
	cusparseDestroy(handle);

	return QMCKL_SUCCESS_DEVICE;
}
#endif

/* mo_select */

// Forward declare this, as its needed by select_mo
qmckl_exit_code_device
qmckl_finalize_mo_basis_device(qmckl_context_device context);

bool qmckl_mo_basis_select_mo_device(qmckl_context_device context,
									 int32_t *keep, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_get_mo_basis_select_mo_device",
									 NULL);
	}

	// WARNING Here, we are expecting a CPU array (instead of a GPU array
	// usually), because it will not be used as a data to be stored in the
	// context. Thus, it makes more sense (and is actually more efficient) to
	// use a CPU array.

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!(qmckl_mo_basis_provided_device(context))) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_mo_basis_select_mo_device",
									 NULL);
	}

	if (keep == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_mo_basis_select_mo_device",
									 "NULL pointer");
	}

	const int64_t mo_num = ctx->mo_basis.mo_num;
	const int64_t ao_num = ctx->ao_basis.ao_num;

	if (size_max < mo_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_mo_basis_select_mo",
									 "Array too small: expected mo_num.");
	}

	int64_t mo_num_new = 0;
	for (int64_t i = 0; i < mo_num; ++i) {
		if (keep[i] != 0)
			++mo_num_new;
	}

	double *restrict coefficient = (double *)qmckl_malloc_device(
		context, ao_num * mo_num_new * sizeof(double));

	int64_t k = 0;
	for (int64_t i = 0; i < mo_num; ++i) {
		if (keep[i] != 0) {
			qmckl_memcpy_D2D(context, &(coefficient[k * ao_num]),
							 &(ctx->mo_basis.coefficient[i * ao_num]),
							 ao_num * sizeof(double));
			++k;
		}
	}

	qmckl_exit_code_device rc =
		qmckl_free_device(context, ctx->mo_basis.coefficient);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->mo_basis.coefficient = coefficient;
	ctx->mo_basis.mo_num = mo_num_new;

	rc = qmckl_finalize_mo_basis_device(context);
	return rc;
}

//**********
// PROVIDE
//**********

/* mo_vgl */

qmckl_exit_code_device
qmckl_provide_mo_basis_mo_vgl_device(qmckl_context_device context) {

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_provide_mo_basis_mo_vgl_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->mo_basis.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_mo_basis_mo_vgl_device",
									 NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->mo_basis.mo_vgl_date) {

		/* Allocate array */
		if (ctx->mo_basis.mo_vgl == NULL) {

			double *mo_vgl = (double *)qmckl_malloc_device(
				context,
				5 * ctx->mo_basis.mo_num * ctx->point.num * sizeof(double));

			if (mo_vgl == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "qmckl_mo_basis_mo_vgl", NULL);
			}
			ctx->mo_basis.mo_vgl = mo_vgl;
		}

		rc = qmckl_provide_ao_basis_ao_vgl_device(context);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
										 "qmckl_ao_basis", NULL);
		}

#if HAVE_CUSPARSE
		rc = qmckl_compute_mo_basis_mo_vgl_cusparse_device(
			context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num, ctx->point.num,
			ctx->mo_basis.coefficient_t, ctx->ao_basis.ao_vgl,
			ctx->mo_basis.mo_vgl);
#elif HAVE_CUBLAS
#if HAVE_FLOATMOS
		rc = qmckl_compute_mo_basis_mo_vgl_sgemm_device(
			context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num, ctx->point.num,
			ctx->mo_basis.coefficient_t, ctx->ao_basis.ao_vgl,
			ctx->mo_basis.mo_vgl);
#else
		rc = qmckl_compute_mo_basis_mo_vgl_dgemm_device(
			context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num, ctx->point.num,
			ctx->mo_basis.coefficient_t, ctx->ao_basis.ao_vgl,
			ctx->mo_basis.mo_vgl);
#endif
#else
		rc = qmckl_compute_mo_basis_mo_vgl_device(
			context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num, ctx->point.num,
			ctx->mo_basis.coefficient_t, ctx->ao_basis.ao_vgl,
			ctx->mo_basis.mo_vgl);
#endif

		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->mo_basis.mo_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

/* mo_value */

qmckl_exit_code_device
qmckl_provide_mo_basis_mo_value_device(qmckl_context_device context) {

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_provide_mo_basis_mo_value_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->mo_basis.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_mo_basis_mo_value_device",
									 NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->mo_basis.mo_value_date) {

		qmckl_exit_code_device rc;

		/* Allocate array */
		if (ctx->mo_basis.mo_value == NULL) {

			double *mo_value = (double *)qmckl_malloc_device(
				context,
				ctx->mo_basis.mo_num * ctx->point.num * sizeof(double));

			if (mo_value == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "qmckl_mo_basis_mo_value", NULL);
			}
			ctx->mo_basis.mo_value = mo_value;
		}

		if (ctx->mo_basis.mo_vgl_date == ctx->point.date) {

			// mo_vgl has been computed at this step: Just copy the data.

			double *v = &(ctx->mo_basis.mo_value[0]);
			double *vgl = &(ctx->mo_basis.mo_vgl[0]);
#pragma omp target is_device_ptr(v, vgl)
			{
				for (int i = 0; i < ctx->point.num; ++i) {
					for (int k = 0; k < ctx->mo_basis.mo_num; ++k) {
						v[k] = vgl[k];
					}
					v += ctx->mo_basis.mo_num;
					vgl += ctx->mo_basis.mo_num * 5;
				}
			}

		} else {

			rc = qmckl_provide_ao_basis_ao_value_device(context);
			if (rc != QMCKL_SUCCESS_DEVICE) {
				return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
											 "qmckl_ao_basis_ao_value_device",
											 NULL);
			}
#if HAVE_CUSPARSE
			rc = qmckl_compute_mo_basis_mo_value_cusparse_device(
				context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num,
				ctx->point.num, ctx->mo_basis.coefficient_t,
				ctx->ao_basis.ao_value, ctx->mo_basis.mo_value);
#elif HAVE_CUBLAS
#if HAVE_FLOATMOS
			rc = qmckl_compute_mo_basis_mo_value_sgemm_device(
				context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num,
				ctx->point.num, ctx->mo_basis.coefficient_t,
				ctx->ao_basis.ao_value, ctx->mo_basis.mo_value);
#else
			rc = qmckl_compute_mo_basis_mo_value_dgemm_device(
				context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num,
				ctx->point.num, ctx->mo_basis.coefficient_t,
				ctx->ao_basis.ao_value, ctx->mo_basis.mo_value);
#endif
#else
			rc = qmckl_compute_mo_basis_mo_value_device(
				context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num,
				ctx->point.num, ctx->mo_basis.coefficient_t,
				ctx->ao_basis.ao_value, ctx->mo_basis.mo_value);
#endif
		}

		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->mo_basis.mo_value_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GET
//**********

/* mo_vgl */

qmckl_exit_code_device
qmckl_get_mo_basis_mo_vgl_device(qmckl_context_device context,
								 double *const mo_vgl, const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_mo_basis_mo_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int64_t sze = 5 * ctx->point.num * ctx->mo_basis.mo_num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_mo_basis_mo_vgl",
									 "input array too small");
	}
	qmckl_memcpy_D2D(context, mo_vgl, ctx->mo_basis.mo_vgl,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_mo_basis_mo_vgl_inplace_device(qmckl_context_device context,
										 double *const mo_vgl,
										 const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device((qmckl_context_device)context,
									 QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_get_mo_basis_mo_vgl_device", NULL);
	}

	qmckl_exit_code_device rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	const int64_t sze = 5 * ctx->mo_basis.mo_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_mo_basis_mo_vgl_device",
									 "input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	double *old_array = ctx->mo_basis.mo_vgl;

	ctx->mo_basis.mo_vgl = mo_vgl;

	rc = qmckl_provide_mo_basis_mo_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->mo_basis.mo_vgl = old_array;

	return QMCKL_SUCCESS_DEVICE;
}

/* mo_value */

qmckl_exit_code_device
qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
								   double *const mo_value,
								   const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_mo_basis_mo_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->point.num * ctx->mo_basis.mo_num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_mo_basis_mo_value",
									 "input array too small");
	}
	qmckl_memcpy_D2D(context, mo_value, ctx->mo_basis.mo_value,
					 sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *const mo_value,
										   const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_NULL_CONTEXT_DEVICE,
			"qmckl_get_mo_basis_mo_value_device", NULL);
	}

	qmckl_exit_code_device rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->mo_basis.mo_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_mo_basis_mo_value_device",
									 "input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	double *old_array = ctx->mo_basis.mo_value;

	ctx->mo_basis.mo_value = mo_value;

	rc = qmckl_provide_mo_basis_mo_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->mo_basis.mo_value = old_array;

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// VARIOUS GETTERS/SETTERS
//**********

qmckl_exit_code_device
qmckl_get_mo_basis_mo_num_device(const qmckl_context_device context,
								 int64_t *mo_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_mo_basis_mo_num", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->mo_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_mo_basis_mo_num", NULL);
	}

	assert(ctx->mo_basis.mo_num > (int64_t)0);
	*mo_num = ctx->mo_basis.mo_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_mo_basis_mo_num_device(qmckl_context_device context, int64_t mo_num) {

	int32_t mask = 1;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_mo_basis_mo_num_device", NULL);
	}

	if (mo_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_mo_basis_mo_num_device",
									 "mo_num <= 0");
	}

	ctx->mo_basis.mo_num = mo_num;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_mo_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_NULL_CONTEXT_DEVICE;
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_mo_basis_coefficient_device",
									 NULL);
	}

	if (ctx->mo_basis.coefficient != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->mo_basis.coefficient);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_mo_basis_coefficient_device", NULL);
		}
	}

	double *new_array = (double *)qmckl_malloc_device(
		context, ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double));
	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_mo_basis_coefficient_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, coefficient,
					 ctx->ao_basis.ao_num * ctx->mo_basis.mo_num *
						 sizeof(double));

	ctx->mo_basis.coefficient = new_array;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}
