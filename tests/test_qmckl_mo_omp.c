#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include "chbrclf.h"
#include "../include/qmckl_gpu.h"
#include <omp.h>

#define MO_VALUE_ID(x, y) chbrclf_mo_num *x + y
#define MO_VGL_ID(x, y, z) 5 * chbrclf_mo_num *x + chbrclf_mo_num *y + z

int main() {
	qmckl_context_device context;

	if (omp_get_num_devices() <= 0) {
		printf("Error : no device found. Aborting execution\n");
		exit(1);
	}
	context = qmckl_context_create_device(omp_get_default_device());

#define walk_num chbrclf_walk_num
#define elec_num chbrclf_elec_num
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num
#define prim_num chbrclf_prim_num

	// Put nucleus stuff in CPU arrays
	int64_t elec_up_num = chbrclf_elec_up_num;
	int64_t elec_dn_num = chbrclf_elec_dn_num;
	double *elec_coord = &(chbrclf_elec_coord[0][0][0]);
	int64_t nucl_num = chbrclf_nucl_num;
	double *nucl_charge = chbrclf_charge;
	double *nucl_coord = &(chbrclf_nucl_coord[0][0]);
	int64_t point_num = walk_num * elec_num;

	// Put nucleus stuff in GPU arrays
	double *elec_coord_d =
		qmckl_malloc_device(context, point_num * 3 * sizeof(double));
	double *nucl_charge_d =
		qmckl_malloc_device(context, nucl_num * sizeof(double));
	double *nucl_coord_d =
		qmckl_malloc_device(context, nucl_num * 3 * sizeof(double));

	// Set nucleus stuff in context

	qmckl_exit_code rc;

	rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);
	assert(rc == QMCKL_SUCCESS);

	assert(qmckl_electron_provided(context));

	rc = qmckl_set_point_device(context, 'N', point_num, elec_coord_d,
								point_num * 3);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_set_nucleus_num_device(context, nucl_num);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord_d,
										nucl_num * 3);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_set_nucleus_charge_device(context, nucl_charge_d, nucl_num);
	assert(rc == QMCKL_SUCCESS);

	assert(qmckl_nucleus_provided(context));

	// Put other stuff in CPU arrays
	int64_t *nucleus_index = &(chbrclf_basis_nucleus_index[0]);
	int64_t *nucleus_shell_num = &(chbrclf_basis_nucleus_shell_num[0]);
	int32_t *shell_ang_mom = &(chbrclf_basis_shell_ang_mom[0]);
	int64_t *shell_prim_num = &(chbrclf_basis_shell_prim_num[0]);
	int64_t *shell_prim_index = &(chbrclf_basis_shell_prim_index[0]);
	double *shell_factor = &(chbrclf_basis_shell_factor[0]);
	double *exponent = &(chbrclf_basis_exponent[0]);
	double *coefficient = &(chbrclf_basis_coefficient[0]);
	double *prim_factor = &(chbrclf_basis_prim_factor[0]);
	double *ao_factor = &(chbrclf_basis_ao_factor[0]);

	// Put other stuff in GPU arrays
	int64_t *nucleus_index_d =
		qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	int64_t *nucleus_shell_num_d =
		qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	int32_t *shell_ang_mom_d =
		qmckl_malloc_device(context, shell_num * sizeof(int32_t));
	int64_t *shell_prim_num_d =
		qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	int64_t *shell_prim_index_d =
		qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	double *shell_factor_d =
		qmckl_malloc_device(context, shell_num * sizeof(double));
	double *exponent_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	double *coefficient_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	double *prim_factor_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	double *ao_factor_d =
		qmckl_malloc_device(context, ao_num * sizeof(double));

	qmckl_memcpy_H2D(context, nucleus_index_d, nucleus_index,
					 nucl_num * sizeof(int64_t));
	qmckl_memcpy_H2D(context, nucleus_shell_num_d, nucleus_shell_num,
					 nucl_num * sizeof(int64_t));
	qmckl_memcpy_H2D(context, shell_ang_mom_d, shell_ang_mom,
					 shell_num * sizeof(int32_t));
	qmckl_memcpy_H2D(context, shell_prim_num_d, shell_prim_num,
					 shell_num * sizeof(int64_t));
	qmckl_memcpy_H2D(context, shell_prim_index_d, shell_prim_index,
					 shell_num * sizeof(int64_t));
	qmckl_memcpy_H2D(context, shell_factor_d, shell_factor,
					 shell_num * sizeof(double));
	qmckl_memcpy_H2D(context, exponent_d, exponent, prim_num * sizeof(double));
	qmckl_memcpy_H2D(context, coefficient_d, coefficient,
					 prim_num * sizeof(double));
	qmckl_memcpy_H2D(context, prim_factor_d, prim_factor,
					 prim_num * sizeof(double));
	qmckl_memcpy_H2D(context, ao_factor_d, ao_factor, ao_num * sizeof(double));

	char typ = 'G';

	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_type_device(context, typ);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_shell_num_device(context, chbrclf_shell_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_prim_num_device(context, chbrclf_prim_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
												 nucl_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_nucleus_shell_num_device(
		context, nucleus_shell_num_d, nucl_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom_d,
												 chbrclf_shell_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor_d,
												chbrclf_shell_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num_d,
												  chbrclf_shell_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index_d,
													chbrclf_shell_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_exponent_device(context, exponent_d,
											chbrclf_prim_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_coefficient_device(context, coefficient_d,
											   chbrclf_prim_num);
	assert(rc == QMCKL_SUCCESS);
	assert(!qmckl_ao_basis_provided(context));

	rc = qmckl_set_ao_basis_prim_factor_device(context, prim_factor_d,
											   chbrclf_prim_num);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_set_ao_basis_ao_num(context, chbrclf_ao_num);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_set_ao_basis_ao_factor_device(context, ao_factor_d,
											 chbrclf_ao_num);
	assert(rc == QMCKL_SUCCESS);

	assert(qmckl_ao_basis_provided(context));

	double *ao_vgl =
		qmckl_malloc_host(context, point_num * 5 * ao_num * sizeof(double));
	double *ao_vgl_d =
		qmckl_malloc_device(context, point_num * 5 * ao_num * sizeof(double));

	double *ao_value =
		qmckl_malloc_host(context, point_num * 5 * ao_num * sizeof(double));
	double *ao_value_d =
		qmckl_malloc_device(context, point_num * 5 * ao_num * sizeof(double));

	rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_d,
										  (int64_t)5 * point_num * ao_num);
	qmckl_memcpy_D2H(context, ao_vgl, ao_vgl_d,
					 point_num * 5 * ao_num * sizeof(double));


	// We thus make sure that ao_value is set
	rc = qmckl_get_ao_basis_ao_value_device(context, ao_vgl_d,
										  (int64_t)point_num * ao_num);
	qmckl_memcpy_D2H(context, ao_value, ao_value_d,
					 point_num * ao_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	/* Set up MO data */
	int64_t mo_num = chbrclf_mo_num;
	rc = qmckl_set_mo_basis_mo_num(context, mo_num);
	assert(rc == QMCKL_SUCCESS);

	double *mo_coefficient = &(chbrclf_mo_coef[0]);
	double *mo_coefficient_d =
		qmckl_malloc_device(context, mo_num * ao_num * sizeof(double));

	rc = qmckl_set_mo_basis_coefficient_device(context, mo_coefficient_d);
	assert(rc == QMCKL_SUCCESS);

	assert(qmckl_mo_basis_provided(context));

	rc = qmckl_context_touch(context);
	assert(rc == QMCKL_SUCCESS);

	double *mo_value =
		qmckl_malloc_host(context, point_num * chbrclf_mo_num * sizeof(double));
	double *mo_value_d = qmckl_malloc_device(
		context, point_num * chbrclf_mo_num * sizeof(double));
	rc = qmckl_get_mo_basis_mo_value_device(context, mo_value_d,
											point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	double * mo_vgl = qmckl_malloc_host(context, point_num * 5 * chbrclf_mo_num);
	double * mo_vgl_d = qmckl_malloc_device(context, point_num * 5 * chbrclf_mo_num);
	rc = qmckl_get_mo_basis_mo_vgl_device(context, &(mo_vgl_d[0]),
								   point_num * 5 * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_vgl, mo_vgl_d,
					 point_num * 5 * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			assert(fabs(mo_vgl[MO_VGL_ID(i, 0, k)] - mo_value[MO_VALUE_ID(i, k)]) < 1.e-12);
		}
	}

	rc = qmckl_context_touch(context);
	assert(rc == QMCKL_SUCCESS);

	rc = qmckl_get_mo_basis_mo_value_device (context, mo_value_d,
									 point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			if(fabs(mo_vgl[MO_VGL_ID(i, 0, k)] - mo_value[MO_VALUE_ID(i, k)]) > 1.e-12) {
				// break and return
			}
		}
	}

	rc = qmckl_get_mo_basis_mo_value_device(context, mo_value,
									 point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);
	#pragma omp target is_device_ptr(mo_vgl, mo_value)
	{
	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			if(fabs(2.0 * mo_vgl[MO_VGL_ID(i, 0, k)] - mo_value[MO_VALUE_ID(i, k)]) >= 1.e-12) {
				// break and return
			};
		}
	}
	}

	// TODO ?
	// rc = qmckl_mo_basis_rescale(context, 0.5);
	assert(rc == QMCKL_SUCCESS);


	printf("\n");
	printf(" mo_vgl mo_vgl[0][26][219] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 219)]);
	printf(" mo_vgl mo_vgl[1][26][219] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 219)]);
	printf(" mo_vgl mo_vgl[0][26][220] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 220)]);
	printf(" mo_vgl mo_vgl[1][26][220] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 220)]);
	printf(" mo_vgl mo_vgl[0][26][221] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 221)]);
	printf(" mo_vgl mo_vgl[1][26][221] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 221)]);
	printf(" mo_vgl mo_vgl[0][26][222] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 222)]);
	printf(" mo_vgl mo_vgl[1][26][222] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 222)]);
	printf(" mo_vgl mo_vgl[0][26][223] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 223)]);
	printf(" mo_vgl mo_vgl[1][26][223] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 223)]);
	printf(" mo_vgl mo_vgl[0][26][224] %25.15e\n", mo_vgl[MO_VGL_ID(0, 26, 224)]);
	printf(" mo_vgl mo_vgl[1][26][224] %25.15e\n", mo_vgl[MO_VGL_ID(1, 26, 224)]);

	// printf(" mo_vgl mo_vgl[0][26][224] %25.15e\n", mo_vgl[2][0][3]);
	// printf(" mo_vgl mo_vgl[1][26][224] %25.15e\n", mo_vgl[2][1][3]);

	printf("\n");

	// TODO
	// rc = qmckl_get_mo_basis_mo_num_device(context, &mo_num);
	// printf(" mo_num: %ld\n", mo_num);
	// assert(mo_num == 2);


	rc = qmckl_context_destroy_device(context);
	assert(rc == QMCKL_SUCCESS);

	return 0;
}
