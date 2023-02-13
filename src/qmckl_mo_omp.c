#include "include/qmckl_mo.h"

//**********
// COMPUTE
//**********

/* mo_vgl */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_device(
	qmckl_context context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t,
	double *restrict ao_vgl, double *restrict mo_vgl) {
	assert(context != QMCKL_NULL_CONTEXT);

#pragma omp target is_device_ptr(coefficient_t, ao_vgl, mo_vgl)
	{
		#pragma omp teams distribute parallel for simd
		for (int64_t j = 0; j < point_num; ++j) {

			// Set j subarray to 0
			for(int k=0; k<5; ++k) {
				for(int l=0; l<mo_num; l++) {
					mo_vgl[l + mo_num*k + mo_num*5*j] = 0.;
				}
			}

			for(int64_t k=0; k<ao_num; k++) {
				if(ao_vgl[k + ao_num*5*j] != 0.) {
					double c1 = ao_vgl[k + ao_num*0 + ao_num*5*j];
					double c2 = ao_vgl[k + ao_num*1 + ao_num*5*j];
					double c3 = ao_vgl[k + ao_num*2 + ao_num*5*j];
					double c4 = ao_vgl[k + ao_num*3 + ao_num*5*j];
					double c5 = ao_vgl[k + ao_num*4 + ao_num*5*j];

					for(int i=0; i<mo_num; i++) {
						mo_vgl[i + mo_num*0 + mo_num*5*j] =
							mo_vgl[i + mo_num*0 + mo_num*5*j] + coefficient_t[i + mo_num*k] * c1;
						mo_vgl[i + mo_num*1 + mo_num*5*j] =
							mo_vgl[i + mo_num*1 + mo_num*5*j] + coefficient_t[i + mo_num*k] * c2;
						mo_vgl[i + mo_num*2 + mo_num*5*j] =
							mo_vgl[i + mo_num*2 + mo_num*5*j] + coefficient_t[i + mo_num*k] * c3;
						mo_vgl[i + mo_num*3 + mo_num*5*j] =
							mo_vgl[i + mo_num*3 + mo_num*5*j] + coefficient_t[i + mo_num*k] * c4;
						mo_vgl[i + mo_num*4 + mo_num*5*j] =
							mo_vgl[i + mo_num*4 + mo_num*5*j] + coefficient_t[i + mo_num*k] * c5;
					}
				}
			}
		}
	}
	// End of GPU region

	return QMCKL_SUCCESS;
}


/* mo_value */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_device(
	qmckl_context context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t,
	double *restrict ao_value, double *restrict mo_value) {
	assert(context != QMCKL_NULL_CONTEXT);

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

	mem_info.size = point_num * ao_num * sizeof(double);
	double *av1_shared = qmckl_malloc_device(context, mem_info);

	mem_info.size = point_num * ao_num * sizeof(int64_t);
	int64_t *idx_shared = qmckl_malloc_device(context, mem_info);

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
				double *restrict ck2 =
					coefficient_t + idx[n + 1] * mo_num;
				double *restrict ck3 =
					coefficient_t + idx[n + 2] * mo_num;
				double *restrict ck4 =
					coefficient_t + idx[n + 3] * mo_num;

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
	return QMCKL_SUCCESS;
}
