#include <qmckl.h>

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_blas_private_type.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"

#include "qmckl_blas_private_func.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_device_types.h"
#include "qmckl_memory_device.h"

qmckl_vector qmckl_vector_alloc_omp_device(qmckl_context_device context,
                                       const int64_t size, int device_id);

qmckl_exit_code qmckl_vector_free_omp_device(qmckl_context_device context,
                                         qmckl_vector *vector, int device_id);

qmckl_exit_code
qmckl_vector_of_double_omp_device(const qmckl_context_device context,
                                  const double *target, const int64_t size_max,
                                  qmckl_vector *vector_out, int device_id);

qmckl_matrix qmckl_matrix_alloc_omp_device(qmckl_context_device context,
                                           const int64_t size1, const int64_t size2,
                                           int device_id);

qmckl_exit_code qmckl_matrix_free_omp_device(qmckl_context_device context,
                                             qmckl_matrix *matrix, int device_id);

qmckl_matrix qmckl_matrix_set_omp_device(qmckl_matrix matrix, double value);

qmckl_exit_code
qmckl_matrix_of_double_omp_device(const qmckl_context_device context,
                                  const double *target, const int64_t size_max,
                                  qmckl_matrix *matrix_out, int device_id);

qmckl_exit_code qmckl_transpose_omp_device(qmckl_context_device context,
                                           const qmckl_matrix A, qmckl_matrix At);

qmckl_tensor qmckl_tensor_alloc_omp_device(qmckl_context context,
                                           const int64_t order, const int64_t *size,
                                           int device_id);

qmckl_exit_code qmckl_tensor_free_omp_device(qmckl_context_device context,
                                             qmckl_tensor *tensor, int device_id);

qmckl_tensor qmckl_tensor_set_omp_device(qmckl_tensor tensor, double value);
