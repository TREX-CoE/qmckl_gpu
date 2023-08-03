#include <stdlib.h>
#include <stdbool.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_memory.h"
qmckl_exit_code_device qmckl_woodbury_kxk(
	const qmckl_context_device context, cublasHandle_t b_handle,
	cusolverDnHandle_t s_handle, const uint64_t Lds, const uint64_t Dim,
	const uint64_t N_updates, const double *__restrict Updates,
	const uint64_t *__restrict Updates_index, const double breakdown,
	double *__restrict Slater_inv, double *__restrict determinant);
