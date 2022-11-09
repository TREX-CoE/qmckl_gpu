#ifndef QMCKL_MO_HPT
#define QMCKL_MO_HPT

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

#include <omp.h>

#include "qmckl_blas_private_type.h"

/* Data structure */

typedef struct qmckl_mo_basis_struct {
	int64_t mo_num;
	double *restrict coefficient;
	double *restrict coefficient_t;

	double *restrict mo_vgl;
	double *restrict mo_value;
	uint64_t mo_vgl_date;
	uint64_t mo_value_date;

	int32_t uninitialized;
	bool provided;

#ifdef HAVE_DEVICE_POINTERS
	double *restrict coefficient_device;
	double *restrict coefficient_t_device;

	double *restrict mo_vgl_device;
	double *restrict mo_value_device;
#endif
} qmckl_mo_basis_struct;

#endif
