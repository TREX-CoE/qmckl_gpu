#ifndef QMCKL_BLAS_HPF
#define QMCKL_BLAS_HPF
#include "qmckl_blas_private_type.h"

qmckl_vector qmckl_vector_alloc(qmckl_context context, const int64_t size);

qmckl_vector qmckl_vector_alloc_device(qmckl_context context,
				       const int64_t size, int device_id);

qmckl_exit_code qmckl_vector_free(qmckl_context context, qmckl_vector *vector);

qmckl_exit_code qmckl_vector_free_device(qmckl_context context,
					 qmckl_vector *vector, int device_id);

qmckl_matrix qmckl_matrix_alloc(qmckl_context context, const int64_t size1,
				const int64_t size2);

#ifdef HAVE_DEVICE_POINTERS
qmckl_matrix qmckl_matrix_alloc_device(qmckl_context context,
				       const int64_t size1, const int64_t size2,
				       int device_id);
#endif

qmckl_exit_code qmckl_matrix_free(qmckl_context context, qmckl_matrix *matrix);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_matrix_free_device(qmckl_context context,
					 qmckl_matrix *matrix, int device_id);
#endif

qmckl_tensor qmckl_tensor_alloc(qmckl_context context, const int64_t order,
				const int64_t *size);

#ifdef HAVE_DEVICE_POINTERS
qmckl_tensor qmckl_tensor_alloc_device(qmckl_context context,
				       const int64_t order, const int64_t *size,
				       int device_id);
#endif

qmckl_exit_code qmckl_tensor_free(qmckl_context context, qmckl_tensor *tensor);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_tensor_free_device(qmckl_context context,
					 qmckl_tensor *tensor, int device_id);
#endif

/* Vector -> Matrix */

qmckl_matrix qmckl_matrix_of_vector(const qmckl_vector vector,
				    const int64_t size1, const int64_t size2);

/* Vector -> Tensor */

qmckl_tensor qmckl_tensor_of_vector(const qmckl_vector vector,
				    const int64_t order, const int64_t *size);

/* Matrix -> Vector */

qmckl_vector qmckl_vector_of_matrix(const qmckl_matrix matrix);

/* Matrix -> Tensor */

qmckl_tensor qmckl_tensor_of_matrix(const qmckl_matrix matrix,
				    const int64_t order, const int64_t *size);

/* Tensor -> Vector */

qmckl_vector qmckl_vector_of_tensor(const qmckl_tensor tensor);

/* Tensor -> Matrix */

qmckl_matrix qmckl_matrix_of_tensor(const qmckl_tensor tensor,
				    const int64_t size1, const int64_t size2);

#define qmckl_vec(v, i) v.data[i]
#define qmckl_mat(m, i, j) m.data[(i) + (j)*m.size[0]]

#define qmckl_ten3(t, i, j, k) t.data[(i) + t.size[0] * ((j) + t.size[1] * (k))]
#define qmckl_ten4(t, i, j, k, l)                                              \
	t.data[(i) + t.size[0] * ((j) + t.size[1] * ((k) + t.size[2] * (l)))]
#define qmckl_ten5(t, i, j, k, l, m)                                           \
	t.data[(i) +                                                           \
	       t.size[0] *                                                     \
		   ((j) +                                                      \
		    t.size[1] * ((k) + t.size[2] * ((l) + t.size[3] * (m))))]

/* Vector */

qmckl_vector qmckl_vector_set(qmckl_vector vector, double value);

/* Matrix */

qmckl_matrix qmckl_matrix_set(qmckl_matrix matrix, double value);

#ifdef HAVE_DEVICE_POINTERS
qmckl_matrix qmckl_matrix_set_device(qmckl_matrix matrix, double value);
#endif

/* Tensor */

qmckl_tensor qmckl_tensor_set(qmckl_tensor tensor, double value);

#ifdef HAVE_DEVICE_POINTERS
qmckl_tensor qmckl_tensor_set_device(qmckl_tensor tensor, double value);
#endif

/* Copy to/from to ~double*~ */

qmckl_exit_code qmckl_double_of_vector(const qmckl_context context,
				       const qmckl_vector vector,
				       double *const target,
				       const int64_t size_max);

qmckl_exit_code qmckl_double_of_matrix(const qmckl_context context,
				       const qmckl_matrix matrix,
				       double *const target,
				       const int64_t size_max);

qmckl_exit_code qmckl_double_of_tensor(const qmckl_context context,
				       const qmckl_tensor tensor,
				       double *const target,
				       const int64_t size_max);

qmckl_exit_code qmckl_vector_of_double(const qmckl_context context,
				       const double *target,
				       const int64_t size_max,
				       qmckl_vector *vector);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_vector_of_double_device(const qmckl_context context,
					      const double *target,
					      const int64_t size_max,
					      qmckl_vector *vector,
					      int device_id);
#endif

qmckl_exit_code qmckl_matrix_of_double(const qmckl_context context,
				       const double *target,
				       const int64_t size_max,
				       qmckl_matrix *matrix);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_matrix_of_double_device(const qmckl_context context,
					      const double *target,
					      const int64_t size_max,
					      qmckl_matrix *matrix,
					      int device_id);
#endif

qmckl_exit_code qmckl_tensor_of_double(const qmckl_context context,
				       const double *target,
				       const int64_t size_max,
				       qmckl_tensor *tensor);

/* ~qmckl_matmul~ */

/*    Matrix multiplication using the =qmckl_matrix= data type: */

/*    \[ */
/*    C_{ij} = \beta C_{ij} + \alpha \sum_{k} A_{ik} \cdot B_{kj} */
/*    \] */

/* #  TODO: Add description about the external library dependence. */

/*    #+NAME: qmckl_matmul_args */
/*    | Variable  | Type            | In/Out | Description       | */
/*    |-----------+-----------------+--------+-------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state      | */
/*    | ~TransA~  | ~char~          | in     | 'T' is transposed | */
/*    | ~TransB~  | ~char~          | in     | 'T' is transposed | */
/*    | ~alpha~   | ~double~        | in     | \alpha            | */
/*    | ~A~       | ~qmckl_matrix~  | in     | Matrix $A$        | */
/*    | ~B~       | ~qmckl_matrix~  | in     | Matrix $B$        | */
/*    | ~beta~    | ~double~        | in     | \beta             | */
/*    | ~C~       | ~qmckl_matrix~  | out    | Matrix $C$        | */

/*    #+CALL:
 * generate_c_header(table=qmckl_matmul_args,rettyp="qmckl_exit_code",fname="qmckl_matmul")
 */

/*     #+RESULTS: */

qmckl_exit_code qmckl_matmul(const qmckl_context context, const char TransA,
			     const char TransB, const double alpha,
			     const qmckl_matrix A, const qmckl_matrix B,
			     const double beta, qmckl_matrix *const C);

/* ~qmckl_transpose~ */

/*    Transposes a matrix: $A^\dagger_{ji} = A_{ij}$. */

/*    | Variable  | Type            | In/Out | Description       | */
/*    |-----------+-----------------+--------+-------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state      | */
/*    | ~A~       | ~qmckl_matrix~  | in     | Input matrix      | */
/*    | ~At~      | ~qmckl_matrix~  | out    | Transposed matrix | */

qmckl_exit_code qmckl_transpose(qmckl_context context, const qmckl_matrix A,
				qmckl_matrix At);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_transpose_device(qmckl_context context,
				       const qmckl_matrix A, qmckl_matrix At);
#endif

#endif
