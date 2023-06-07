#include "../include/qmckl_distance.h"

qmckl_exit_code_device
qmckl_distance_device(const qmckl_context_device context, const char transa,
					  const char transb, const int64_t m, const int64_t n,
					  const double *A, const int64_t lda, const double *B,
					  const int64_t ldb, double *const C, const int64_t ldc) {
	int i, j;
	double x, y, z;
	int transab;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (m <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (n <= 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

	if (transa == 'N' || transa == 'n') {
		transab = 0;
	} else if (transa == 'T' || transa == 't') {
		transab = 1;
	} else {
		transab = -100;
	}

	if (transb == 'N' || transb == 'n') {
	} else if (transb == 'T' || transb == 't') {
		transab = transab + 2;
	} else {
		transab = -100;
	}

	if (transab < 0) {
		info = QMCKL_INVALID_ARG_1_DEVICE;
		return info;
	}

	// check for LDA
	if ((transab & 1) == 0 && lda < 3) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 1) == 1 && lda < m) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2) == 0 && lda < 3) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2) == 2 && lda < m) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	// check for LDB
	if ((transab & 1) == 0 && ldb < 3) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 1) == 1 && ldb < n) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 0 && ldb < 3) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 2 && ldb < n) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	// check for LDC
	if (ldc < m) {
		info = QMCKL_INVALID_ARG_11_DEVICE;
		return info;
	}

	switch (transab) {

	case 0:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[0 + j * ldb];
				y = A[1 + i * lda] - B[1 + j * ldb];
				z = A[2 + i * lda] - B[2 + j * ldb];
				C[i + j * ldc] = x * x + y * y + z * z;
			}
			for (int i = 0; i < ldc; i++)
				C[i + j * ldc] = sqrt(C[i + j * ldc]);
		}
	} break;

	case 1:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[0 + j * ldb];
				y = A[i + 1 * lda] - B[1 + j * ldb];
				z = A[i + 2 * lda] - B[2 + j * ldb];
				C[i + j * ldc] = x * x + y * y + z * z;
			}
			for (int i = 0; i < j; i++)
				C[i + j * ldc] = sqrt(C[i + j * ldc]);
		}
	} break;

	case 2:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[j + 0 * ldb];
				y = A[1 + i * lda] - B[j + 1 * ldb];
				z = A[2 + i * lda] - B[j + 2 * ldb];
				C[i + j * ldc] = x * x + y * y + z * z;
			}
			for (int i = 0; i < ldc; i++)
				C[i + j * ldc] = sqrt(C[i + j * ldc]);
		}
	} break;

	case 3:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[j + 0 * ldb];
				y = A[i + 1 * lda] - B[j + 1 * ldb];
				z = A[i + 2 * lda] - B[j + 2 * ldb];

				C[i + j * ldc] = x * x + y * y + z * z;
			}
			for (int i = 0; i < ldc; i++)
				C[i + j * ldc] = sqrt(C[i + j * ldc]);
		}
	} break;
	}

	return info;
}

qmckl_exit_code_device qmckl_distance_rescaled_device(
	const qmckl_context_device context, const char transa, const char transb,
	const int64_t m, const int64_t n, const double *A, const int64_t lda,
	const double *B, const int64_t ldb, double *const C, const int64_t ldc,
	const double rescale_factor_kappa) {

	int i, j, transab;
	double x, y, z, dist, rescale_factor_kappa_inv;

	rescale_factor_kappa_inv = 1.0 / rescale_factor_kappa;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (m <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (n <= 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

	if (transa == 'N' || transa == 'n') {
		transab = 0;
	} else if (transa == 'T' || transa == 't') {
		transab = 1;
	} else {
		transab = -100;
	}

	if (transb == 'N' || transb == 'n') {
	} else if (transb == 'T' || transb == 't') {
		transab = transab + 2;
	} else {
		transab = -100;
	}

	// check for LDA
	if (transab < 0) {
		info = QMCKL_INVALID_ARG_1_DEVICE;
		return info;
	}

	if ((transab & 1) == 0 && lda < 3) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 1) == 1 && lda < m) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2) == 0 && lda < 3) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2) == 2 && lda < m) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	// check for LDB
	if ((transab & 1) == 0 && ldb < 3) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 1) == 1 && ldb < n) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 0 && ldb < 3) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 2 && ldb < n) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	// check for LDC
	if (ldc < m) {
		info = QMCKL_INVALID_ARG_11_DEVICE;
		return info;
	}

	switch (transab) {

	case 0:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[0 + j * ldb];
				y = A[1 + i * lda] - B[1 + j * ldb];
				z = A[2 + i * lda] - B[2 + j * ldb];
				dist = sqrt(x * x + y * y + z * z);
				C[i + j * ldc] = (1.0 - exp(-rescale_factor_kappa * dist)) *
								 rescale_factor_kappa_inv;
			}
		}
	} break;

	case 1:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[0 + j * ldb];
				y = A[i + 1 * lda] - B[1 + j * ldb];
				z = A[i + 2 * lda] - B[2 + j * ldb];
				dist = sqrt(x * x + y * y + z * z);
				C[i + j * ldc] = (1.0 - exp(-rescale_factor_kappa * dist)) *
								 rescale_factor_kappa_inv;
			}
		}
	} break;

	case 2:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[j + 0 * ldb];
				y = A[1 + i * lda] - B[j + 1 * ldb];
				z = A[2 + i * lda] - B[j + 2 * ldb];
				dist = sqrt(x * x + y * y + z * z);
				C[i + j * ldc] = (1.0 - exp(-rescale_factor_kappa * dist)) *
								 rescale_factor_kappa_inv;
			}
		}
	} break;

	case 3:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[j + 0 * ldb];
				y = A[i + 1 * lda] - B[j + 1 * ldb];
				z = A[i + 2 * lda] - B[j + 2 * ldb];
				dist = sqrt(x * x + y * y + z * z);
				C[i + j * ldc] = (1.0 - exp(-rescale_factor_kappa * dist)) *
								 rescale_factor_kappa_inv;
			}
		}
	} break;
	}

	return info;
}

qmckl_exit_code_device qmckl_distance_rescaled_deriv_e_device(
	qmckl_context_device context, char transa, char transb, int m, int n,
	double *A, int lda, double *B, int ldb, double *C, int ldc,
	double rescale_factor_kappa) {

	double x, y, z, dist, dist_inv;
	double rescale_factor_kappa_inv, rij;
	int transab;

	rescale_factor_kappa_inv = 1.0 / rescale_factor_kappa;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (m <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (n <= 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

	if (transa == 'N' || transa == 'n') {
		transab = 0;
	} else if (transa == 'T' || transa == 't') {
		transab = 1;
	} else {
		transab = -100;
	}

	if (transb == 'N' || transb == 'n') {
	} else if (transb == 'T' || transb == 't') {
		transab = transab + 2;
	} else {
		transab = -100;
	}

	// check for LDA
	if (transab < 0) {
		info = QMCKL_INVALID_ARG_1_DEVICE;
		return info;
	}

	if ((transab & 1 == 0) && (lda < 3)) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 1 == 1) && (lda < m)) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2 == 0) && (lda < 3)) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	if ((transab & 2 == 2) && (lda < m)) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	// check for LDB
	if ((transab & 1) == 0 && (ldb < 3)) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 1) == 1 && (ldb < n)) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 0 && (ldb < 3)) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	if ((transab & 2) == 2 && (ldb < n)) {
		info = QMCKL_INVALID_ARG_9_DEVICE;
		return info;
	}

	// check for LDC
	if (ldc < m) {
		info = QMCKL_INVALID_ARG_11_DEVICE;
		return info;
	}

	switch (transab) {

	case 0:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[0 + j * ldb];
				y = A[1 + i * lda] - B[1 + j * ldb];
				z = A[2 + i * lda] - B[2 + j * ldb];
				dist = sqrt(x * x + y * y + z * z);
				dist_inv = 1.0 / dist;
				rij = (1.0 - exp(-rescale_factor_kappa * dist)) *
					  rescale_factor_kappa_inv;
				C[0 + i * 4 + j * 4 * ldc] =
					x * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[1 + i * 4 + j * 4 * ldc] =
					y * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[2 + i * 4 + j * 4 * ldc] =
					z * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[3 + i * 4 + j * 4 * ldc] =
					(2.0 * dist_inv - rescale_factor_kappa_inv) *
					(1.0 - rescale_factor_kappa_inv * rij);
			}
		}
	} break;

	case 1:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[0 + j * ldb];
				y = A[i + 1 * lda] - B[1 + j * ldb];
				z = A[i + 2 * lda] - B[2 + j * ldb];
				dist = sqrt(x * x + y * y + z * z);
				dist_inv = 1.0 / dist;
				rij = (1.0 - exp(-rescale_factor_kappa * dist)) *
					  rescale_factor_kappa_inv;
				C[0 + i + j] =
					x * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[1 + i + j] =
					y * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[2 + i + j] =
					z * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[3 + i + j] = (2.0 * dist_inv - rescale_factor_kappa_inv) *
							   (1.0 - rescale_factor_kappa_inv * rij);
			}
		}
	} break;

	case 2:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[0 + i * lda] - B[j + 0 * ldb];
				y = A[1 + i * lda] - B[j + 1 * ldb];
				z = A[2 + i * lda] - B[j + 2 * ldb];
				dist = sqrt(x * x + y * y + z * z);
				dist_inv = 1.0 / dist;
				rij = (1.0 - exp(-rescale_factor_kappa * dist)) *
					  rescale_factor_kappa_inv;
				C[0 + i * 4 + j * 4 * ldc] =
					x * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[1 + i * 4 + j * 4 * ldc] =
					y * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[2 + i * 4 + j * 4 * ldc] =
					z * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[3 + i * 4 + j * ldc] =
					(2.0 * dist_inv - rescale_factor_kappa_inv) *
					(1.0 - rescale_factor_kappa_inv * rij);
			}
		}
	} break;

	case 3:
#pragma omp target is_device_ptr(A, B, C)
	{
#pragma omp teams distribute parallel for simd collapse(2)
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				x = A[i + 0 * lda] - B[j + 0 * ldb];
				y = A[i + 1 * lda] - B[j + 1 * ldb];
				z = A[i + 2 * lda] - B[j + 2 * ldb];
				dist = sqrt(x * x + y * y + z * z);
				dist_inv = 1.0 / dist;
				rij = (1.0 - exp(-rescale_factor_kappa * dist)) *
					  rescale_factor_kappa_inv;
				C[0 + i * 4 + j * 4 * ldc] =
					x * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[1 + i * 4 + j * 4 * ldc] =
					y * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[2 + i * 4 + j * 4 * ldc] =
					z * dist_inv * (1.0 - rescale_factor_kappa_inv * rij);
				C[3 + i * 4 + j * 4 * ldc] =
					(2.0 * dist_inv - rescale_factor_kappa_inv) *
					(1.0 - rescale_factor_kappa_inv * rij);
			}
		}
	} break;
	}

	return info;
}
