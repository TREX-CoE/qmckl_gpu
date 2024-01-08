#include "../include/qmckl_mo.h"

//**********
// COMPUTE
//**********

/* mo_vgl */

qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t, double *restrict ao_vgl,
	double *restrict mo_vgl) {

	assert(context != QMCKL_NULL_CONTEXT_DEVICE);


#ifdef HAVE_LIBGPUBLAS

		gpu_dgemm('N', 'N', mo_num, point_num*5, ao_num, 1., coefficient_t, mo_num, ao_vgl, ao_num, 0., mo_vgl, mo_num);

#else

#pragma omp target is_device_ptr(coefficient_t, ao_vgl, mo_vgl)
	{

#pragma omp teams distribute parallel for
		for (int64_t j = 0; j < point_num; ++j) {

			// Set j subarray to 0
#pragma omp simd
			for (int k = 0; k < 5; ++k) {
				for (int l = 0; l < mo_num; l++) {
					mo_vgl[l + mo_num * k + mo_num * 5 * j] = 0.;
				}
			}

			for (int64_t k = 0; k < ao_num; k++) {
				for (int l = 0; l < 5; l++) {
					if (ao_vgl[k + ao_num * 5 * j] != 0.) {
						//double c1 = ao_vgl[k + ao_num * l + ao_num * 5 * j];
						double c1 = ao_vgl[k + ao_num * l + ao_num * 5 * j];
					
#pragma omp simd

						for (int i = 0; i < mo_num; i++) {
				/*			mo_vgl[i + mo_num * l + mo_num * 5 * j] =
								mo_vgl[i + mo_num * l + mo_num * 5 * j] +
								coefficient_t[i + mo_num * k] * c1; */
							mo_vgl[j + point_num * l + point_num * 5 * i] =
								mo_vgl[j + point_num * l + point_num * 5 * i] +
								coefficient_t[k + ao_num * i] * c1;
				

						}
					}
				}
			}
		}
	}
#endif

	// End of GPU region

	return QMCKL_SUCCESS_DEVICE;
}

/* mo_value */

qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t,
	double *restrict ao_value, double *restrict mo_value) {
	assert(context != QMCKL_NULL_CONTEXT_DEVICE);

	double *av1_shared =
		qmckl_malloc_device(context, point_num * ao_num * sizeof(double));
	int64_t *idx_shared =
		qmckl_malloc_device(context, point_num * ao_num * sizeof(int64_t));

#pragma omp target is_device_ptr(coefficient_t, ao_value, mo_value,            \
								 idx_shared, av1_shared)
	{
#pragma omp teams distribute parallel for
		for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {

			double *av1 = av1_shared + ipoint * ao_num;
			int64_t *idx = idx_shared + ipoint * ao_num;

			double *vgl1 = mo_value + ipoint * mo_num;
			double *avgl1 = ao_value + ipoint * ao_num;

			for (int64_t i = 0; i < mo_num; ++i) {
				vgl1[i] = 0.;
			}

			int64_t nidx = 0;
			for (int64_t k = 0; k < ao_num; ++k) {
				if (avgl1[k] != 0.) {
					idx[nidx] = k;
					av1[nidx] = avgl1[k];
					++nidx;
				}
			}

			int64_t n = 0;

			for (n = 0; n < nidx - 4; n += 4) {
				double *restrict ck1 = coefficient_t + idx[n] * mo_num;
				double *restrict ck2 = coefficient_t + idx[n + 1] * mo_num;
				double *restrict ck3 = coefficient_t + idx[n + 2] * mo_num;
				double *restrict ck4 = coefficient_t + idx[n + 3] * mo_num;

				double a11 = av1[n];
				double a21 = av1[n + 1];
				double a31 = av1[n + 2];
				double a41 = av1[n + 3];

#pragma omp simd
				for (int64_t i = 0; i < mo_num; ++i) {
					vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 +
							  ck3[i] * a31 + ck4[i] * a41;
				}
			}

			for (int64_t m = n; m < nidx; m += 1) {
				double *restrict ck = coefficient_t + idx[m] * mo_num;
				double a1 = av1[m];

#pragma omp simd
				for (int64_t i = 0; i < mo_num; ++i) {
					vgl1[i] += ck[i] * a1;
				}
			}
		}
	}
	return QMCKL_SUCCESS_DEVICE;
}

//**********
// FINALIZE MO BASIS
//**********

qmckl_exit_code_device
qmckl_finalize_mo_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_finalize_mo_basis_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	double *new_array = (double *)qmckl_malloc_device(
		context, ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double));
	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_finalize_mo_basis_device", NULL);
	}

	assert(ctx->mo_basis.coefficient != NULL);

	if (ctx->mo_basis.coefficient_t != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->mo_basis.coefficient_t);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_finalize_mo_basis_device", NULL);
		}
	}

	double *coefficient = ctx->mo_basis.coefficient;

	int64_t ao_num = ctx->ao_basis.ao_num;
	int64_t mo_num = ctx->mo_basis.mo_num;

#pragma omp target is_device_ptr(new_array, coefficient)
	{
#pragma omp parallel for collapse(2)
		for (int64_t i = 0; i < ao_num; ++i) {
			for (int64_t j = 0; j < mo_num; ++j) {
				new_array[i * mo_num + j] = coefficient[j * ao_num + i];
			}
		}
	}

	ctx->mo_basis.coefficient_t = new_array;
	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;
	return rc;
}
