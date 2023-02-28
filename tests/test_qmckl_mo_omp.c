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

#define walk_num chbrclf_walk_num
#define elec_num chbrclf_elec_num
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num
#define prim_num chbrclf_prim_num

#define AO_VALUE_ID(x, y) ao_num *x + y
#define MO_VALUE_ID(x, y) chbrclf_mo_num *x + y
#define MO_VGL_ID(x, y, z) 5 * chbrclf_mo_num *x + chbrclf_mo_num *y + z
#define AO_VGL_ID(x, y, z) 5 * ao_num *x + ao_num *y + z

int main() {
	qmckl_context_device context;

	if (omp_get_num_devices() <= 0) {
		printf("Error : no device found. Aborting execution\n");
		exit(1);
	}
	context = qmckl_context_create_device(omp_get_default_device());

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

	qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge,
					 nucl_num * sizeof(double));
	qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord,
					 3 * nucl_num * sizeof(double));
	qmckl_memcpy_H2D(context, elec_coord_d, elec_coord,
					 3 * point_num * sizeof(double));

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
	double *ao_factor_d = qmckl_malloc_device(context, ao_num * sizeof(double));

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

	double *ao_vgl = malloc(point_num * 5 * ao_num * sizeof(double));
	double *ao_vgl_d =
		qmckl_malloc_device(context, point_num * 5 * ao_num * sizeof(double));

	double *ao_value = malloc(point_num * ao_num * sizeof(double));
	double *ao_value_d =
		qmckl_malloc_device(context, point_num * ao_num * sizeof(double));

	rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_d,
										  (int64_t)5 * point_num * ao_num);
	qmckl_memcpy_D2H(context, ao_vgl, ao_vgl_d,
					 point_num * 5 * ao_num * sizeof(double));

	/* From here, context is completely initialized, we just
	 * make sure that the ao_value array is set (because we
	 * have computed ao_vgl) */

	rc = qmckl_get_ao_basis_ao_value_device(context, ao_value_d,
											(int64_t)point_num * ao_num);
	qmckl_memcpy_D2H(context, ao_value, ao_value_d,
					 point_num * ao_num * sizeof(double));

	/* Check values from AO, so as to make sure that if we have errors in
	 * MO later, they do come from MO */

	printf("\n");
	printf(" ao_vgl ao_vgl[26][0][219] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 219)]);
	printf(" ao_vgl ao_vgl[26][1][219] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 219)]);
	printf(" ao_vgl ao_vgl[26][2][219] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 219)]);
	printf(" ao_vgl ao_vgl[26][3][219] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 219)]);
	printf(" ao_vgl ao_vgl[26][4][219] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 219)]);
	printf("\n");

	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 219)] - (1.020298798341620e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 219)] - (-4.928035238010602e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 219)] - (-4.691009312035986e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 219)] - (1.449504046436699e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 219)] - (4.296442111843973e-07)) > 1.e-14)
		return 1;

	printf(" ao_value ao_value[26][219] %25.15e\n",
		   ao_value[AO_VALUE_ID(26, 219)]);
	printf(" ao_value ao_value[26][220] %25.15e\n",
		   ao_value[AO_VALUE_ID(26, 220)]);
	printf(" ao_value ao_value[26][221] %25.15e\n",
		   ao_value[AO_VALUE_ID(26, 221)]);
	printf(" ao_value ao_value[26][222] %25.15e\n",
		   ao_value[AO_VALUE_ID(26, 222)]);
	printf("\n");

	if (fabs(ao_value[AO_VALUE_ID(26, 219)] - (1.020298798341620e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_value[AO_VALUE_ID(26, 220)] - (1.516643537739178e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_value[AO_VALUE_ID(26, 221)] - (-4.686370882518819e-09)) >
		1.e-14)
		return 1;
	if (fabs(ao_value[AO_VALUE_ID(26, 222)] - (7.514816980753531e-09)) > 1.e-14)
		return 1;

	assert(rc == QMCKL_SUCCESS);

	/* Set up MO data */

	int64_t mo_num = chbrclf_mo_num;
	rc = qmckl_set_mo_basis_mo_num(context, mo_num);
	assert(rc == QMCKL_SUCCESS);

	double *mo_coefficient = &(chbrclf_mo_coef[0]);
	double *mo_coefficient_d =
		qmckl_malloc_device(context, mo_num * ao_num * sizeof(double));
	qmckl_memcpy_H2D(context, mo_coefficient_d, mo_coefficient,
					 mo_num * ao_num * sizeof(double));

	rc = qmckl_set_mo_basis_coefficient_device(context, mo_coefficient_d);
	assert(rc == QMCKL_SUCCESS);

	assert(qmckl_mo_basis_provided(context));

	rc = qmckl_context_touch(context);
	assert(rc == QMCKL_SUCCESS);

	/* Get MO value (from scratch) */
	double *mo_value = malloc(point_num * chbrclf_mo_num * sizeof(double));
	double *mo_value_d = qmckl_malloc_device(
		context, point_num * chbrclf_mo_num * sizeof(double));
	rc = qmckl_get_mo_basis_mo_value_device(context, mo_value_d,
											point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	/* Get MO vgl */
	double *mo_vgl = malloc(point_num * 5 * chbrclf_mo_num * sizeof(double));
	double *mo_vgl_d = qmckl_malloc_device(
		context, point_num * 5 * chbrclf_mo_num * sizeof(double));
	rc = qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl_d,
										  point_num * 5 * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_vgl, mo_vgl_d,
					 point_num * 5 * chbrclf_mo_num * sizeof(double));

	// Making sure that value element of vgl == value
	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			if (fabs(mo_vgl[MO_VGL_ID(i, 0, k)] - mo_value[MO_VALUE_ID(i, k)]) >
				1.e-12) {
				return 1;
			}
		}
	}

	rc = qmckl_context_touch(context);
	assert(rc == QMCKL_SUCCESS);

	/* Get MO value (from MO vgl array) */
	rc = qmckl_get_mo_basis_mo_value_device(context, mo_value_d,
											point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	// Making sure that value element of vgl == value
	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			if (fabs(mo_vgl[MO_VGL_ID(i, 0, k)] - mo_value[MO_VALUE_ID(i, k)]) >
				1.e-12) {
				return 1;
			}
		}
	}

	rc = qmckl_get_mo_basis_mo_value_device(context, mo_value_d,
											point_num * chbrclf_mo_num);
	qmckl_memcpy_D2H(context, mo_value, mo_value_d,
					 point_num * chbrclf_mo_num * sizeof(double));

	assert(rc == QMCKL_SUCCESS);

	// Making sure that value element of vgl == value
	for (int i = 0; i < point_num; ++i) {
		for (int k = 0; k < chbrclf_mo_num; ++k) {
			if (fabs(mo_vgl[MO_VGL_ID(i, 0, k)] -
					 mo_value[MO_VALUE_ID(i, k)]) >= 1.e-12) {
				return 1;
			};
		}
	}

	// TODO Add support for mo_basis_rescale at some point ?
	// rc = qmckl_mo_basis_rescale(context, 0.5);
	// assert(rc == QMCKL_SUCCESS);

	printf("\n");
	printf(" mo_vgl mo_vgl[0][0][0] %25.15e\n", mo_vgl[MO_VGL_ID(0, 0, 0)]);
	printf(" mo_vgl mo_vgl[0][0][1] %25.15e\n", mo_vgl[MO_VGL_ID(0, 0, 1)]);
	printf(" mo_vgl mo_vgl[0][0][2] %25.15e\n", mo_vgl[MO_VGL_ID(0, 0, 2)]);
	printf(" mo_vgl mo_vgl[0][0][3] %25.15e\n", mo_vgl[MO_VGL_ID(0, 0, 3)]);
	printf(" mo_vgl mo_vgl[0][1][0] %25.15e\n", mo_vgl[MO_VGL_ID(0, 1, 0)]);
	printf(" mo_vgl mo_vgl[0][1][1] %25.15e\n", mo_vgl[MO_VGL_ID(0, 1, 1)]);
	printf(" mo_vgl mo_vgl[0][1][2] %25.15e\n", mo_vgl[MO_VGL_ID(0, 1, 2)]);
	printf(" mo_vgl mo_vgl[0][1][3] %25.15e\n", mo_vgl[MO_VGL_ID(0, 1, 3)]);
	printf(" mo_vgl mo_vgl[0][4][0] %25.15e\n", mo_vgl[MO_VGL_ID(0, 4, 0)]);
	printf(" mo_vgl mo_vgl[0][4][1] %25.15e\n", mo_vgl[MO_VGL_ID(0, 4, 1)]);
	printf(" mo_vgl mo_vgl[0][4][2] %25.15e\n", mo_vgl[MO_VGL_ID(0, 4, 2)]);
	printf(" mo_vgl mo_vgl[0][4][3] %25.15e\n", mo_vgl[MO_VGL_ID(0, 4, 3)]);
	printf("\n");


	// Read the ao_vgl ref

	// We will try to open ao_reference.txt "from" qmckl_gpu/ and qmckl_gpu/tests/
	FILE* fp = fopen("tests/mo_reference.txt", "r");
	if (fp == NULL) {
		fp = fopen("mo_reference.txt", "r");
	}
	if(fp == NULL) {
		printf("Error : mo_reference.txt not found, leaving\n");
		exit(1);
	}

	double ref;
	for(int i=0; i<point_num; i++) {
		for(int j=0; j<5; j++) {
			for(int k=0; k<mo_num; k++) {
				fscanf(fp, "%lf", &ref);
				if (fabs(mo_vgl[MO_VGL_ID(i, j, k)] - ref) > 1.e-14)
					return 1;
			}
		}
	}


	if (fabs(mo_vgl[MO_VGL_ID(0, 4, 3)] - (1.848343348520137e-06)) > 1.e-14)
		return 1;

	// TODO
	// rc = qmckl_get_mo_basis_mo_num_device(context, &mo_num);
	// printf(" mo_num: %ld\n", mo_num);
	// assert(mo_num == 2);

	// TODO
	// rc = qmckl_context_destroy_device(context);
	assert(rc == QMCKL_SUCCESS);

	return 0;
}
