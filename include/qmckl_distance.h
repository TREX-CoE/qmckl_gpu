#include <math.h>

#include "qmckl_types.h"

qmckl_exit_code_device qmckl_distance_rescaled_device(
	const qmckl_context_device context, const char transa, const char transb,
	const int64_t m, const int64_t n, const double *A, const int64_t lda,
	const double *B, const int64_t ldb, double *const C, const int64_t ldc,
	const double rescale_factor_kappa);
