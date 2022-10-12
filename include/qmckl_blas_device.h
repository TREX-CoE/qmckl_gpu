#include <qmckl.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_blas_private_type.h"

#include "qmckl_memory_private_func.h"
#include "qmckl_blas_private_func.h"

#include "qmckl_device_types.h"
#include "qmckl_memory_device.h"


qmckl_vector_device
qmckl_vector_alloc_device( qmckl_context_device context,
                           const int64_t size,
                           int device_id);

qmckl_exit_code
qmckl_vector_free_device( qmckl_context_device context,
                          qmckl_vector_device* vector,
                          int device_id);

qmckl_exit_code
qmckl_vector_of_double_device(const qmckl_context_device context,
                              const double* target,
                              const int64_t size_max,
                              qmckl_vector_device* vector_out,
                              int device_id);

qmckl_matrix_device
qmckl_matrix_alloc_device( qmckl_context_device context,
                           const int64_t size1,
                           const int64_t size2,
                           int device_id);

qmckl_exit_code
qmckl_matrix_free_device( qmckl_context_device context,
                          qmckl_matrix_device* matrix,
                          int device_id);

qmckl_matrix_device
qmckl_matrix_set_device(qmckl_matrix_device matrix, double value);

qmckl_exit_code
qmckl_matrix_of_double_device(const qmckl_context_device context,
                              const double* target,
                              const int64_t size_max,
                              qmckl_matrix_device* matrix_out,
                              int device_id);

qmckl_exit_code
qmckl_transpose_device (qmckl_context_device context,
                        const qmckl_matrix_device A,
                        qmckl_matrix_device At );

qmckl_tensor_device
qmckl_tensor_alloc_device( qmckl_context_device context,
                           const int64_t  order,
                           const int64_t* size,
                           int device_id);

qmckl_exit_code
qmckl_tensor_free_device( qmckl_context_device context,
                          qmckl_tensor_device* tensor,
                          int device_id);

qmckl_tensor_device
qmckl_tensor_set_device(qmckl_tensor_device tensor, double value);
