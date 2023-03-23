#include "include/qmckl_mo.h"

#ifdef HAVE_CUBLAS 
#include <cublas_v2.h>
#endif

#ifdef HAVE_NVBLAS 
#include <nvblas.h>
#endif

#ifdef HAVE_CUSPARSE 
#include <cuda_runtime_api.h>
#include <cusparse_v2.h>
#endif

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_acc_offload(qmckl_context context,
													  double *const mo_vgl,
													  const int64_t size_max) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_mo_basis_mo_vgl_acc_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->point.num * 5 * ctx->mo_basis.mo_num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_vgl_acc_offload",
							  "input array too small");
	}
	memcpy(mo_vgl, ctx->mo_basis.mo_vgl, sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_mo_basis_mo_vgl_acc_offload_inplace(
	qmckl_context context, double *const mo_vgl, const int64_t size_max) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_mo_basis_mo_vgl_acc_offload", NULL);
	}

	qmckl_exit_code rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	const int64_t sze = ctx->mo_basis.mo_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_mo_basis_mo_vgl_acc_offload",
							  "input array too small");
	}

	rc = qmckl_context_touch(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	double *old_array = ctx->mo_basis.mo_vgl;

	ctx->mo_basis.mo_vgl = mo_vgl;

	rc = qmckl_provide_mo_basis_mo_vgl_acc_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->mo_basis.mo_vgl = old_array;

	return QMCKL_SUCCESS;
}

//======================================================================
qmckl_exit_code
qmckl_provide_mo_basis_mo_vgl_acc_offload(qmckl_context context) {

	qmckl_exit_code rc = QMCKL_SUCCESS;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_mo_basis_mo_vgl_acc_offload",
							  NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->mo_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_provide_mo_basis_mo_vgl_acc_offload",
							  NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->mo_basis.mo_vgl_date) {

		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size =
			5 * ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

		if (ctx->mo_basis.mo_vgl != NULL) {
			qmckl_memory_info_struct mem_info_test =
				qmckl_memory_info_struct_zero;
			rc = qmckl_get_malloc_info(context, ctx->mo_basis.mo_vgl,
									   &mem_info_test);

			/* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function
			   because the memory was not allocated with qmckl_malloc */

			if ((rc == QMCKL_SUCCESS) &&
				(mem_info_test.size != mem_info.size)) {
				rc = qmckl_free(context, ctx->mo_basis.mo_vgl);
				assert(rc == QMCKL_SUCCESS);
				ctx->mo_basis.mo_vgl = NULL;
			}
		}

		/* Allocate array */
		if (ctx->mo_basis.mo_vgl == NULL) {

			double *mo_vgl = (double *)qmckl_malloc(context, mem_info);

			if (mo_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_mo_basis_mo_vgl", NULL);
			}
			ctx->mo_basis.mo_vgl = mo_vgl;
		}

		rc = qmckl_provide_ao_vgl_acc_offload(context);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith(context, QMCKL_NOT_PROVIDED, "qmckl_ao_basis",
								  NULL);
		}

#pragma acc enter data copyin(ctx)
#pragma acc enter data copyin(ctx->ao_basis.ao_num)
#pragma acc enter data copyin(ctx->mo_basis.mo_num)
#pragma acc enter data copyin(ctx->point.num)
#pragma acc enter data copyin(                                                 \
		ctx->mo_basis                                                          \
			.coefficient_t[0 : ctx->ao_basis.ao_num * ctx->mo_basis.mo_num])
#pragma acc enter data create(                                                 \
		ctx->mo_basis.mo_vgl[0 : ctx->point.num * 5 * ctx->mo_basis.mo_num])

		rc = qmckl_compute_mo_basis_mo_vgl_acc_offload(
			context, ctx->ao_basis.ao_num, ctx->mo_basis.mo_num, ctx->point.num,
			ctx->mo_basis.coefficient_t, ctx->ao_basis.ao_vgl,
			ctx->mo_basis.mo_vgl);

#pragma acc exit data copyout(                                                 \
		ctx->mo_basis.mo_vgl[0 : ctx->point.num * 5 * ctx->mo_basis.mo_num])
#pragma acc exit data copyout(ctx)

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->mo_basis.mo_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *coefficient_t, const double *ao_vgl,
	double *const mo_vgl) {
#if defined(HAVE_CUBLAS)
	return qmckl_compute_mo_basis_mo_vgl_cublas_acc_offload(
		context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#elif defined(HAVE_CUSPARSE)
	return qmckl_compute_mo_basis_mo_vgl_cusparse_acc_offload(
		context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#elif defined(HAVE_NVBLAS)
	return qmckl_compute_mo_basis_mo_vgl_nvblas_acc_offload(
		context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#else
	return qmckl_compute_mo_basis_mo_vgl_hpc_acc_offload(
		context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
	// return qmckl_compute_mo_basis_mo_vgl_doc_acc_offload (context, ao_num,
	// mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#endif
}

#if defined(HAVE_CUSPARSE)
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_cusparse_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

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

#pragma acc host_data use_device(ao_vgl[point_num * 5 * ao_num],               \
									 mo_vgl[point_num * 5 * mo_num],           \
									 coefficient_t[ao_num * mo_num])
	{
		cusparseCreateDnMat(&matAd, A_nrows, A_ncols, lda, ao_vgl,
							cuda_datatype, A_storage_order);
		cusparseCreateDnMat(&matBd, B_nrows, B_ncols, ldb, coefficient_t,
							cuda_datatype, B_storage_order);
		cusparseCreateDnMat(&matCd, C_nrows, C_ncols, ldc, mo_vgl,
							cuda_datatype, C_storage_order);
	}

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
	// int64_t matsize = 5*point_num*ao_num
	// printf("ao_vgl zeros = %d out of %d elements, sparse % = %f\n",
	//       (matsize-As_nnz), matsize, ((float) matsize-As_nnz)/matsize);

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

	cusparseDestroyDnMat(matAd);
	cusparseDestroyDnMat(matBd);
	cusparseDestroyDnMat(matCd);
	cusparseDestroySpMat(matAs);
	cusparseDestroy(handle);

	return QMCKL_SUCCESS;
}
#endif

#if defined(HAVE_NVBLAS)
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_nvblas_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

	double const alpha = 1.;
	double const beta = 0.;

	// NB: cublas views arrays as column-major (ie FORTRAN-like) ordered
	int const m = mo_num;
	int const k = ao_num;
	int const n = point_num * 5;

	int const lda = m;
	int const ldb = k;
	int const ldc = m;

#pragma acc host_data use_device(ao_vgl[point_num * 5 * ao_num],               \
									 mo_vgl[point_num * 5 * mo_num],           \
									 coefficient_t[ao_num * mo_num])
	{
		dgemm("N", "N", &m, &n, &k, &alpha, coefficient_t, &lda, ao_vgl, &ldb,
			  &beta, mo_vgl, &ldc);
	}

	return QMCKL_SUCCESS;
}
#endif

#if defined(HAVE_CUBLAS)
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_cublas_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

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
#pragma acc host_data use_device(ao_vgl[point_num * 5 * ao_num],               \
									 mo_vgl[point_num * 5 * mo_num],           \
									 coefficient_t[ao_num * mo_num])
	{
		cublasDgemm_v2(handle, transa, transb, m, n, k, &alpha, coefficient_t,
					   lda, ao_vgl, ldb, &beta, mo_vgl, ldc);
	}
	cublasDestroy(handle);

	return QMCKL_SUCCESS;
}
#endif

//==================================================================================
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_doc_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

// memset(mo_vgl, 0,  5 * mo_num * point_num * sizeof(double));
#pragma acc parallel loop default(present)
	for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {
		for (int64_t icomp = 0; icomp < 5; ++icomp) {
			for (int64_t i = 0; i < mo_num; ++i) {
				mo_vgl[i + icomp * mo_num + ipoint * mo_num * 5] = 0;
			}
		}
	}

#pragma acc parallel loop default(present)
	for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {
		for (int64_t icomp = 0; icomp < 5; ++icomp) {
			for (int64_t k = 0; k < ao_num; ++k) {
				if (ao_vgl[k + ipoint * ao_num * 5] != 0.) {
					double c = ao_vgl[k + icomp * ao_num + ipoint * ao_num * 5];
					for (int64_t i = 0; i < mo_num; ++i) {
						mo_vgl[i + icomp * mo_num + ipoint * mo_num * 5] =
							mo_vgl[i + icomp * mo_num + ipoint * mo_num * 5] +
							coefficient_t[i + k * ao_num] * c;
					}
				}
			}
		}
	}

	return QMCKL_SUCCESS;
}

//==================================================================================
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_hpc_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

	int64_t idx[ao_num];
	double av1[ao_num];
	double av2[ao_num];
	double av3[ao_num];
	double av4[ao_num];
	double av5[ao_num];

	double vgl1[mo_num];
	double vgl2[mo_num];
	double vgl3[mo_num];
	double vgl4[mo_num];
	double vgl5[mo_num];

#pragma acc enter data create(                                                 \
		idx[0 : ao_num], av1[0 : ao_num], av2[0 : ao_num], av3[0 : ao_num],    \
			av4[0 : ao_num], av5[0 : ao_num], vgl1[0 : mo_num],                \
			vgl2[0 : mo_num], vgl3[0 : mo_num], vgl4[0 : mo_num],              \
			vgl5[0 : mo_num])

#pragma acc parallel num_gangs(point_num)                                      \
	present(coefficient_t[0 : ao_num * mo_num],                                \
				ao_vgl[0 : point_num * 5 * ao_num],                            \
				mo_vgl[0 : point_num * 5 * mo_num])
#pragma acc loop gang private(                                                 \
		idx[0 : ao_num], av1[0 : ao_num], av2[0 : ao_num], av3[0 : ao_num],    \
			av4[0 : ao_num], av5[0 : ao_num], vgl1[0 : mo_num],                \
			vgl2[0 : mo_num], vgl3[0 : mo_num], vgl4[0 : mo_num],              \
			vgl5[0 : mo_num])
	for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {

#pragma acc loop vector
		for (int64_t i = 0; i < mo_num; ++i) {
			vgl1[i] = 0.;
			vgl2[i] = 0.;
			vgl3[i] = 0.;
			vgl4[i] = 0.;
			vgl5[i] = 0.;
		}

		int64_t nidx = 0;
		for (int64_t k = 0; k < ao_num; ++k) {
			if (ao_vgl[k + ipoint * 5 * ao_num] != 0.) {
				idx[nidx] = k;
				av1[nidx] = ao_vgl[k + ipoint * 5 * ao_num];
				av2[nidx] = ao_vgl[k + ao_num + ipoint * 5 * ao_num];
				av3[nidx] = ao_vgl[k + 2 * ao_num + ipoint * 5 * ao_num];
				av4[nidx] = ao_vgl[k + 3 * ao_num + ipoint * 5 * ao_num];
				av5[nidx] = ao_vgl[k + 4 * ao_num + ipoint * 5 * ao_num];
				++nidx;
			}
		}

		int64_t n = 0;

		for (n = 0; n < nidx - 4; n += 4) {
			size_t i1 = idx[n] * mo_num;
			size_t i2 = idx[n + 1] * mo_num;
			size_t i3 = idx[n + 2] * mo_num;
			size_t i4 = idx[n + 3] * mo_num;

			const double a11 = av1[n];
			const double a21 = av1[n + 1];
			const double a31 = av1[n + 2];
			const double a41 = av1[n + 3];

			const double a12 = av2[n];
			const double a22 = av2[n + 1];
			const double a32 = av2[n + 2];
			const double a42 = av2[n + 3];

			const double a13 = av3[n];
			const double a23 = av3[n + 1];
			const double a33 = av3[n + 2];
			const double a43 = av3[n + 3];

			const double a14 = av4[n];
			const double a24 = av4[n + 1];
			const double a34 = av4[n + 2];
			const double a44 = av4[n + 3];

			const double a15 = av5[n];
			const double a25 = av5[n + 1];
			const double a35 = av5[n + 2];
			const double a45 = av5[n + 3];

#pragma acc loop vector
			for (int64_t i = 0; i < mo_num; ++i) {
				// TG: += here results in slowdown
				vgl1[i] = vgl1[i] + coefficient_t[i + i1] * a11 +
						  coefficient_t[i + i2] * a21 +
						  coefficient_t[i + i3] * a31 +
						  coefficient_t[i + i4] * a41;
				vgl2[i] = vgl2[i] + coefficient_t[i + i1] * a12 +
						  coefficient_t[i + i2] * a22 +
						  coefficient_t[i + i3] * a32 +
						  coefficient_t[i + i4] * a42;
				vgl3[i] = vgl3[i] + coefficient_t[i + i1] * a13 +
						  coefficient_t[i + i2] * a23 +
						  coefficient_t[i + i3] * a33 +
						  coefficient_t[i + i4] * a43;
				vgl4[i] = vgl4[i] + coefficient_t[i + i1] * a14 +
						  coefficient_t[i + i2] * a24 +
						  coefficient_t[i + i3] * a34 +
						  coefficient_t[i + i4] * a44;
				vgl5[i] = vgl5[i] + coefficient_t[i + i1] * a15 +
						  coefficient_t[i + i2] * a25 +
						  coefficient_t[i + i3] * a35 +
						  coefficient_t[i + i4] * a45;
			}
		}

		for (int64_t m = n; m < nidx; m += 1) {
			size_t ii = idx[m] * mo_num;
			const double a1 = av1[m];
			const double a2 = av2[m];
			const double a3 = av3[m];
			const double a4 = av4[m];
			const double a5 = av5[m];

			for (int64_t i = 0; i < mo_num; ++i) {
				vgl1[i] += coefficient_t[i + ii] * a1;
				vgl2[i] += coefficient_t[i + ii] * a2;
				vgl3[i] += coefficient_t[i + ii] * a3;
				vgl4[i] += coefficient_t[i + ii] * a4;
				vgl5[i] += coefficient_t[i + ii] * a5;
			}
		}

		for (int64_t i = 0; i < mo_num; ++i) {
			mo_vgl[i + ipoint * mo_num * 5] = vgl1[i];
			mo_vgl[i + mo_num + ipoint * mo_num * 5] = vgl2[i];
			mo_vgl[i + 2 * mo_num + ipoint * mo_num * 5] = vgl3[i];
			mo_vgl[i + 3 * mo_num + ipoint * mo_num * 5] = vgl4[i];
			mo_vgl[i + 4 * mo_num + ipoint * mo_num * 5] = vgl5[i];
		}
	}
#pragma acc exit data delete (                                                 \
		vgl1[0 : mo_num], vgl2[0 : mo_num], vgl3[0 : mo_num],                  \
			vgl4[0 : mo_num], vgl5[0 : mo_num], idx[0 : ao_num],               \
			av1[0 : ao_num], av2[0 : ao_num], av3[0 : ao_num],                 \
			av4[0 : ao_num], av5[0 : ao_num])

	return QMCKL_SUCCESS;
}

//==================================================================================
qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_hpc(
	const qmckl_context context, const int64_t ao_num, const int64_t mo_num,
	const int64_t point_num, const double *restrict coefficient_t,
	const double *restrict ao_vgl, double *restrict const mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

	for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {
		double *restrict const vgl1 = &(mo_vgl[ipoint * 5 * mo_num]);
		double *restrict const vgl2 = vgl1 + mo_num;
		double *restrict const vgl3 = vgl1 + (mo_num << 1);
		double *restrict const vgl4 = vgl1 + (mo_num << 1) + mo_num;
		double *restrict const vgl5 = vgl1 + (mo_num << 2);

		const double *restrict avgl1 = &(ao_vgl[ipoint * 5 * ao_num]);
		const double *restrict avgl2 = avgl1 + ao_num;
		const double *restrict avgl3 = avgl1 + (ao_num << 1);
		const double *restrict avgl4 = avgl1 + (ao_num << 1) + ao_num;
		const double *restrict avgl5 = avgl1 + (ao_num << 2);

		for (int64_t i = 0; i < mo_num; ++i) {
			vgl1[i] = 0.;
			vgl2[i] = 0.;
			vgl3[i] = 0.;
			vgl4[i] = 0.;
			vgl5[i] = 0.;
		}

		int64_t nidx = 0;
		int64_t idx[ao_num];
		double av1[ao_num];
		double av2[ao_num];
		double av3[ao_num];
		double av4[ao_num];
		double av5[ao_num];
		for (int64_t k = 0; k < ao_num; ++k) {
			if (avgl1[k] != 0.) {
				idx[nidx] = k;
				av1[nidx] = avgl1[k];
				av2[nidx] = avgl2[k];
				av3[nidx] = avgl3[k];
				av4[nidx] = avgl4[k];
				av5[nidx] = avgl5[k];
				++nidx;
			}
		}

		int64_t n = 0;

		for (n = 0; n < nidx - 4; n += 4) {
			const double *restrict ck1 = coefficient_t + idx[n] * mo_num;
			const double *restrict ck2 = coefficient_t + idx[n + 1] * mo_num;
			const double *restrict ck3 = coefficient_t + idx[n + 2] * mo_num;
			const double *restrict ck4 = coefficient_t + idx[n + 3] * mo_num;

			const double a11 = av1[n];
			const double a21 = av1[n + 1];
			const double a31 = av1[n + 2];
			const double a41 = av1[n + 3];

			const double a12 = av2[n];
			const double a22 = av2[n + 1];
			const double a32 = av2[n + 2];
			const double a42 = av2[n + 3];

			const double a13 = av3[n];
			const double a23 = av3[n + 1];
			const double a33 = av3[n + 2];
			const double a43 = av3[n + 3];

			const double a14 = av4[n];
			const double a24 = av4[n + 1];
			const double a34 = av4[n + 2];
			const double a44 = av4[n + 3];

			const double a15 = av5[n];
			const double a25 = av5[n + 1];
			const double a35 = av5[n + 2];
			const double a45 = av5[n + 3];

#pragma omp simd
			for (int64_t i = 0; i < mo_num; ++i) {
				vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 +
						  ck4[i] * a41;
				vgl2[i] = vgl2[i] + ck1[i] * a12 + ck2[i] * a22 + ck3[i] * a32 +
						  ck4[i] * a42;
				vgl3[i] = vgl3[i] + ck1[i] * a13 + ck2[i] * a23 + ck3[i] * a33 +
						  ck4[i] * a43;
				vgl4[i] = vgl4[i] + ck1[i] * a14 + ck2[i] * a24 + ck3[i] * a34 +
						  ck4[i] * a44;
				vgl5[i] = vgl5[i] + ck1[i] * a15 + ck2[i] * a25 + ck3[i] * a35 +
						  ck4[i] * a45;
			}
		}

		for (int64_t m = n; m < nidx; m += 1) {
			const double *restrict ck = coefficient_t + idx[m] * mo_num;
			const double a1 = av1[m];
			const double a2 = av2[m];
			const double a3 = av3[m];
			const double a4 = av4[m];
			const double a5 = av5[m];

#pragma omp simd
			for (int64_t i = 0; i < mo_num; ++i) {
				vgl1[i] += ck[i] * a1;
				vgl2[i] += ck[i] * a2;
				vgl3[i] += ck[i] * a3;
				vgl4[i] += ck[i] * a4;
				vgl5[i] += ck[i] * a5;
			}
		}
	}

	return QMCKL_SUCCESS;
}
