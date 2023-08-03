#include "include/qmckl_woodbury.h"

qmckl_exit_code_device qmckl_woodbury_kxk(
	const qmckl_context_device context, cublasHandle_t b_handle,
	cusolverDnHandle_t s_handle, const uint64_t Lds, const uint64_t Dim,
	const uint64_t N_updates, const double *__restrict Updates,
	const uint64_t *__restrict Updates_index, const double breakdown,
	double *__restrict Slater_inv, double *__restrict determinant) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith(context, QMCKL_NULL_CONTEXT_DEVICE, "qmckl_woodbury_kxk",
							  NULL);
	}

	bool swap;
	uint32_t j;
	double alpha, beta, det;
	int *__restrict pivot = malloc(sizeof *pivot * N_updates);
	double *__restrict __attribute__((aligned(8))) C =
		malloc(sizeof *C * Dim * N_updates);
	double *__restrict __attribute__((aligned(8))) B =
		malloc(sizeof *B * N_updates * N_updates);
	double *__restrict __attribute__((aligned(8))) Binv =
		malloc(sizeof *Binv * N_updates * N_updates);
	double *__restrict __attribute__((aligned(8))) D =
		malloc(sizeof *D * N_updates * Lds);
	double *__restrict __attribute__((aligned(8))) T1 =
		malloc(sizeof *T1 * N_updates * Lds);
	double *__restrict __attribute__((aligned(8))) T2 =
		malloc(sizeof *T2 * Dim * Lds);

	int workspace_size = 0, info;
	double *__restrict __attribute__((aligned(8))) workspace = NULL;
	cusolverDnDgetrf_bufferSize(s_handle, N_updates, N_updates, B, N_updates,
								&workspace_size);
	workspace = malloc(sizeof *workspace * workspace_size);

#pragma omp target enter data map(to                                           \
								  : Updates [0:Lds * N_updates],               \
									Updates_index [0:N_updates],               \
									Slater_inv [0:Dim * Lds])

	// Compute C <- S^{-1} U : Dim x K : standard dgemm
	alpha = 1.0f, beta = 0.0f;
#pragma omp target enter data map(alloc : C [0:Dim * N_updates])
#pragma omp target data use_device_ptr(Slater_inv, Updates, C)
	{
		(void)cublasDgemm(b_handle, CUBLAS_OP_T, CUBLAS_OP_N, N_updates, Dim,
						  Lds, &alpha, Updates, Lds, Slater_inv, Lds, &beta, C,
						  N_updates);
	}
#pragma omp target exit data map(delete : Updates [0:Lds * N_updates])

// Construct B <- 1 + V C : K x K
#pragma omp target enter data map(alloc : B [0:N_updates * N_updates])
#pragma omp target teams distribute parallel for collapse(2)                   \
	map(present, alloc                                                         \
		: B, C, Updates_index)
	for (int j = 0; j < N_updates; ++j) {
		for (int i = 0; i < N_updates; ++i) {
			const uint32_t row = Updates_index[i] - 1;
			B[i * N_updates + j] = C[row * N_updates + j] + (i == j);
		}
	}

// Compute det(B) via LU(B)
#pragma omp target enter data map(alloc                                        \
								  : workspace [0:workspace_size],              \
									pivot [0:N_updates], info)
#pragma omp target data use_device_ptr(B, workspace, pivot, info)
	{
		(void)cusolverDnDgetrf(s_handle, N_updates, N_updates, B, N_updates,
							   workspace, pivot,
							   &info); // col-maj enforced, so res. is LU(B)^T
	}
#pragma omp target exit data map(delete : workspace [0:workspace_size])
	swap = false;
	j = 0;
	det = 1.0f;
#pragma omp target teams distribute parallel for reduction(+ : j)              \
	reduction(* : det)
	for (uint32_t i = 0; i < N_updates; i++) {
		swap = (bool)(pivot[i] -
					  (i + 1)); // swap = {0->false/no swap, >0->true/swap}
		j += (uint32_t)swap;	// count # of swaps
		det *= B[i * (N_updates + 1)]; // prod. of diag elm. of B
	}
	if (fabs(det) < breakdown)
		return 1; // check if determinant of B is too close to zero. If so, exit
				  // early.
	if (determinant) { // update det(Slater) if determinant != NULL
		if ((j & 1) != 0)
			det = -det; // multiply det with -1 if # of swaps is odd
		*determinant *= det;
	}

// Compute B^{-1} : initialise as I for solving BX=I
#pragma omp target enter data map(alloc : Binv [0:N_updates * N_updates])
#pragma omp target teams distribute parallel for collapse(2)
	for (int i = 0; i < N_updates; ++i) {
		for (int j = 0; j < N_updates; ++j) {
			Binv[i * N_updates + j] = (i == j);
		}
	}

#pragma omp target data use_device_ptr(B, pivot, Binv, info)
	{
		(void)cusolverDnDgetrs(s_handle, CUBLAS_OP_T, N_updates, N_updates, B,
							   N_updates, pivot, Binv, N_updates,
							   &info); // Needs op(B) = B^T because of line 403
	}
#pragma omp target exit data map(delete                                        \
								 : B [0:N_updates * N_updates],                \
								   pivot [0:N_updates], info)

// Construct D = V S^{-1} : K x LDS
#pragma omp target enter data map(alloc : D [0:N_updates * Lds])
#pragma omp target teams distribute parallel for collapse(2)
	for (uint32_t i = 0; i < N_updates; ++i) {
		for (uint32_t j = 0; j < Lds; ++j) {
			const uint32_t row = Updates_index[i] - 1;
			D[i * Lds + j] = Slater_inv[row * Lds + j];
		}
	}
#pragma omp target exit data map(delete : Updates_index [0:N_updates])

// T1 <- B^{-1} D : KxLDS : standard dgemm
#pragma omp target enter data map(alloc : T1 [0:N_updates * Lds])
#pragma omp target data use_device_ptr(D, Binv, T1)
	{
		(void)cublasDgemm(
			b_handle, CUBLAS_OP_N,
			CUBLAS_OP_T, // REMEMBER THIS IS Binv TRANSPOSED  because of
						 // cusolverDnDgetrs CALL ON l.434 !!!
			Lds, N_updates, N_updates, &alpha, D, Lds, Binv, N_updates, &beta,
			T1, Lds);
	}

#pragma omp target exit data map(delete                                        \
								 : D [0:N_updates * Lds],                      \
								   Binv [0:N_updates * N_updates])

	// Compute S^{-1} <- S^{-1} - C * T1 : Dim x LDS : standard dgemm
	alpha = -1.0f, beta = 1.0f;
#pragma omp target data use_device_ptr(T1, C, Slater_inv)
	{
		(void)cublasDgemm(b_handle, CUBLAS_OP_N, CUBLAS_OP_N, Dim, Lds,
						  N_updates, &alpha, T1, Lds, C, N_updates, &beta,
						  Slater_inv, Lds);
	}

#pragma omp target exit data map(delete                                        \
								 : T1 [0:N_updates * Lds],                     \
								   C [0:Dim * N_updates])
#pragma omp target update from(Slater_inv [0:Dim * Lds])
#pragma omp target exit data map(delete : Slater_inv [0:Dim * Lds])

	free(pivot);
	free(B);
	free(Binv);
	free(C);
	free(D);
	free(T1);

	return QMCKL_SUCCESS_DEVICE;
}
