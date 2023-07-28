
#include <rocblas/rocblas.h>
#include <hip/hip_runtime.h>

//*********
// DGEMM
//*********


qmckl_exit_code_device qmckl_gpu_dgemm(qmckl_context_device context, char transA, char transB, int64_t m, int64_t n, int64_t k, double alpha, double* A, int64_t lda,\
 double* B, int64_t ldb, double beta, double* C, int64_t ldc);