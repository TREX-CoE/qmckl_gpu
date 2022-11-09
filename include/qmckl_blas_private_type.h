#ifndef QMCKL_BLAS_HPT
#define QMCKL_BLAS_HPT

/* Vector */

/*   | Variable | Type      | Description             | */
/*   |----------+-----------+-------------------------| */
/*   | ~size~   | ~int64_t~ | Dimension of the vector | */
/*   | ~data~   | ~double*~ | Elements                | */

typedef struct qmckl_vector {
	double *restrict data;

#ifdef HAVE_DEVICE_POINTERS
	double *restrict data_device;
#endif

	int64_t size;
} qmckl_vector;

/* Matrix */

/*   | Variable | Type         | Description                 | */
/*   |----------+--------------+-----------------------------| */
/*   | ~size~   | ~int64_t[2]~ | Dimension of each component | */
/*   | ~data~   | ~double*~    | Elements                    | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */

typedef struct qmckl_matrix {
	double *restrict data;

#ifdef HAVE_DEVICE_POINTERS
	double *restrict data_device;
#endif

	int64_t size[2];
} qmckl_matrix;

/* Tensor */

/*   | Variable | Type                              | Description | */
/*   |----------+-----------------------------------+-----------------------------|
 */
/*   | ~order~  | ~int64_t~                         | Order of the tensor | */
/*   | ~size~   | ~int64_t[QMCKL_TENSOR_ORDER_MAX]~ | Dimension of each
 * component | */
/*   | ~data~   | ~double*~                         | Elements | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */

#define QMCKL_TENSOR_ORDER_MAX 16

typedef struct qmckl_tensor {
	double *restrict data;
	int64_t order;
	int64_t size[QMCKL_TENSOR_ORDER_MAX];

#ifdef HAVE_DEVICE_POINTERS
	double *restrict data_device;
#endif
} qmckl_tensor;

#endif
