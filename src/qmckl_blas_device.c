#include "../include/qmckl_blas_device.h"
#include <omp.h>

// NOTE Device versions of these functions will be added as they are needed

//**********
// VECTOR
//**********

qmckl_vector_device qmckl_vector_alloc_device(qmckl_context_device context,
                                              const int64_t size,
                                              int device_id) {
  /* Should always be true by contruction */
  assert(size > (int64_t)0);

  qmckl_vector_device result;
  result.size = size;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size * sizeof(double);
  result.data_device = qmckl_malloc_device(context, mem_info, device_id);

  if (result.data_device == NULL) {
    result.size = (int64_t)0;
  }

  return result;
}

qmckl_exit_code qmckl_vector_free_device(qmckl_context_device context,
                                         qmckl_vector_device *vector,
                                         int device_id) {
  /* Always true */
  assert(vector->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, vector->data_device, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  vector->size = (int64_t)0;
  vector->data = NULL;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_vector_of_double_device(const qmckl_context_device context,
                              const double *target, const int64_t size_max,
                              qmckl_vector_device *vector_out, int device_id) {

  // Accepts an host array an copies it in the device section of vector_out
  // (assuming the vector is already allocated)

  qmckl_vector_device vector = *vector_out;
  /* Always true by construction */
  assert(qmckl_context_check_device(context) != QMCKL_NULL_CONTEXT);

  if (vector.size == 0) {
    // This error is thrown
    return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
                          "qmckl_vector_of_double", "Vector not allocated");
  }

  if (vector.size != size_max) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
                          "qmckl_vector_of_double", "Wrong vector size");
  }

  omp_target_memcpy(vector.data_device, target, vector.size * sizeof(double), 0,
                    0, device_id, omp_get_initial_device());

  *vector_out = vector;
  return QMCKL_SUCCESS;
}

//**********
// MATRIX
//**********

qmckl_matrix_device qmckl_matrix_alloc_device(qmckl_context_device context,
                                              const int64_t size1,
                                              const int64_t size2,
                                              int device_id) {
  /* Should always be true by contruction */
  assert(size1 * size2 > (int64_t)0);

  qmckl_matrix_device result;

  result.size[0] = size1;
  result.size[1] = size2;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size1 * size2 * sizeof(double);
  result.data_device =
      (double *)qmckl_malloc_device(context, mem_info, device_id);

  if (result.data_device == NULL) {
    result.size[0] = (int64_t)0;
    result.size[1] = (int64_t)0;
  }

  return result;
}

qmckl_exit_code qmckl_matrix_free_device(qmckl_context_device context,
                                         qmckl_matrix_device *matrix,
                                         int device_id) {
  /* Always true */
  assert(matrix->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, matrix->data_device, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }
  matrix->data_device = NULL;
  matrix->size[0] = (int64_t)0;
  matrix->size[1] = (int64_t)0;

  return QMCKL_SUCCESS;
}

qmckl_matrix_device qmckl_matrix_set_device(qmckl_matrix_device matrix,
                                            double value) {
  // Recompute array size
  int prod_size = matrix.size[0] * matrix.size[1];

  double *data_device = matrix.data_device;
#pragma omp target is_device_ptr(data_device)
  {
    for (int i = 0; i < prod_size; i++) {
      data_device[i] = value;
    }
  }
  return matrix;
}

qmckl_exit_code
qmckl_matrix_of_double_device(const qmckl_context_device context,
                              const double *target, const int64_t size_max,
                              qmckl_matrix_device *matrix_out, int device_id) {

  // Accepts an host array an copies it in the device section of matrix
  // (assuming the matrix is already allocated)

  qmckl_matrix_device matrix = *matrix_out;
  /* Always true by construction */
  assert(qmckl_context_check_device(context) != QMCKL_NULL_CONTEXT);

  if (matrix.size[0] * matrix.size[1] == 0) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
                          "qmckl_matrix_of_double", "Matrix not allocated");
  }

  if (matrix.size[0] * matrix.size[1] != size_max) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
                          "qmckl_matrix_of_double", "Wrong vector size");
  }

  omp_target_memcpy(matrix.data_device, target, size_max * sizeof(double), 0, 0,
                    device_id, omp_get_initial_device());

  *matrix_out = matrix;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_transpose_device(qmckl_context_device context,
                                       const qmckl_matrix_device A,
                                       qmckl_matrix_device At) {
  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (A.size[0] < 1) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_transpose",
                          "Invalid size for A");
  }

  if (At.data_device == NULL) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_transpose",
                          "Output matrix not allocated");
  }

  if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
    return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_transpose",
                          "Invalid size for At");
  }

  double *A_data = A.data_device;
  int A_s0 = A.size[0];

  double *At_data = At.data_device;
  int At_s0 = At.size[0];
  int At_s1 = At.size[1];

#pragma omp target is_device_ptr(A_data, At_data)
  {
#pragma omp parallel for collapse(2)
    for (int64_t j = 0; j < At_s1; ++j)
      for (int64_t i = 0; i < At_s0; ++i)
        At_data[i + j * At_s0] = A_data[j + i * A_s0];
  }

  return QMCKL_SUCCESS;
}

//**********
// TENSOR
//**********

qmckl_tensor_device qmckl_tensor_alloc_device(qmckl_context_device context,
                                              const int64_t order,
                                              const int64_t *size,
                                              int device_id) {
  /* Should always be true by contruction */
  assert(order > 0);
  assert(order <= QMCKL_TENSOR_ORDER_MAX);
  assert(size != NULL);

  qmckl_tensor_device result;
  result.order = order;

  int64_t prod_size = (int64_t)1;
  for (int64_t i = 0; i < order; ++i) {
    assert(size[i] > (int64_t)0);
    result.size[i] = size[i];
    prod_size *= size[i];
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prod_size * sizeof(double);

  result.data_device =
      (double *)qmckl_malloc_device(context, mem_info, device_id);

  if (result.data_device == NULL) {
    memset(&result, 0, sizeof(qmckl_tensor));
  }

  return result;
}

qmckl_exit_code qmckl_tensor_free_device(qmckl_context_device context,
                                         qmckl_tensor_device *tensor,
                                         int device_id) {
  /* Always true */
  assert(tensor->data_device != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, tensor->data_device, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  // TODO Memset to 0
  // memset(tensor, 0, sizeof(qmckl_tensor));

  return QMCKL_SUCCESS;
}

qmckl_tensor_device qmckl_tensor_set_device(qmckl_tensor_device tensor,
                                            double value) {
  // Recompute array size
  int prod_size = 1;

  for (int i = 0; i < tensor.order; i++) {
    prod_size *= tensor.size[i];
  }

  double *data_device = tensor.data_device;
#pragma omp target is_device_ptr(data_device)
  {
    for (int i = 0; i < prod_size; i++) {
      data_device[i] = value;
    }
  }
  return tensor;
}
