#include <stdlib.h>
#include <stdbool.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include <qmckl.h>

qmckl_exit_code qmckl_woodbury_kxk(
    cublasHandle_t b_handle,
    cusolverDnHandle_t s_handle,
    const uint64_t Lds,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* __restrict Updates,
    const uint64_t* __restrict Updates_index,
    const double breakdown,
    double* __restrict Slater_inv,
    double* __restrict determinant);
