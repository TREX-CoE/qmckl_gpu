
#include "../include/qmckl_gpu_blas.hpp"
//*********
// DGEMM
//*********


qmckl_exit_code_device qmckl_gpu_dgemm(qmckl_context_device context, char transA, char transB, int64_t m, int64_t n, int64_t k, double alpha, double* A, int64_t lda,\
 double* B, int64_t ldb, double beta, double* C, int64_t ldc)
{
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