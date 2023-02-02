#include "qmckl.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "chbrclf.h"
//#include "../qmckl/src/qmckl_ao_private_func.h"
#include "../include/qmckl_gpu.h"

#define AO_VGL_ID(x, y, z) 5 * ao_num *x + ao_num *y + z

int main() {
	qmckl_context context;
	// TODO Get device ID according to OpenMP/ACC
	context = qmckl_context_create_device(0);

	const int64_t nucl_num = chbrclf_nucl_num;

	// Put nucleus stuff in CPU arrays
	const double *nucl_charge = chbrclf_charge;
	const double *nucl_coord = &(chbrclf_nucl_coord[0][0]);

	// Put nucleus stuff in GPU arrays
	const double *nucl_charge_d =
		qmckl_malloc_device(context, nucl_num * sizeof(double));
	const double *nucl_coord_d =
		qmckl_malloc_device(context, 3 * nucl_num * sizeof(double));

	qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge,
					 nucl_num * sizeof(double));
	qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord,
					 3 * nucl_num * sizeof(double));

	// Set nucleus stuff in context
	qmckl_exit_code rc;
	rc = qmckl_set_nucleus_num_device(context, nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord_d,
										3 * nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_nucleus_charge_device(context, nucl_charge_d, nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	if (!qmckl_nucleus_provided(context))
		return 1;

	const int64_t shell_num = chbrclf_shell_num;
	const int64_t prim_num = chbrclf_prim_num;
	const int64_t ao_num = chbrclf_ao_num;

	// Put other stuff in CPU arrays
	const int64_t *nucleus_index = &(chbrclf_basis_nucleus_index[0]);
	const int64_t *nucleus_shell_num = &(chbrclf_basis_nucleus_shell_num[0]);
	const int32_t *shell_ang_mom = &(chbrclf_basis_shell_ang_mom[0]);
	const int64_t *shell_prim_num = &(chbrclf_basis_shell_prim_num[0]);
	const int64_t *shell_prim_index = &(chbrclf_basis_shell_prim_index[0]);
	const double *shell_factor = &(chbrclf_basis_shell_factor[0]);
	const double *exponent = &(chbrclf_basis_exponent[0]);
	const double *coefficient = &(chbrclf_basis_coefficient[0]);
	const double *prim_factor = &(chbrclf_basis_prim_factor[0]);
	const double *ao_factor = &(chbrclf_basis_ao_factor[0]);

	// Put other stuff in GPU arrays
	const int64_t *nucleus_index_d =
		qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	const int64_t *nucleus_shell_num_d =
		qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	const int32_t *shell_ang_mom_d =
		qmckl_malloc_device(context, shell_num * sizeof(int32_t));
	const int64_t *shell_prim_num_d =
		qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	const int64_t *shell_prim_index_d =
		qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	const double *shell_factor_d =
		qmckl_malloc_device(context, shell_num * sizeof(double));
	const double *exponent_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	const double *coefficient_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	const double *prim_factor_d =
		qmckl_malloc_device(context, prim_num * sizeof(double));
	const double *ao_factor_d =
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

	const char typ = 'G';

	rc = qmckl_set_ao_basis_type_device(context, typ);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
												 nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_nucleus_shell_num_device(
		context, nucleus_shell_num_d, nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom_d,
												 shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor_d,
												shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num_d,
												  shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index_d,
													shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_exponent_device(context, exponent_d, prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc =
		qmckl_set_ao_basis_coefficient_device(context, coefficient_d, prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc =
		qmckl_set_ao_basis_prim_factor_device(context, prim_factor_d, prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	rc = qmckl_set_ao_basis_ao_factor_device(context, ao_factor_d, ao_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	if (!qmckl_ao_basis_provided(context))
		return 1;

	// Checking arrays after context set and get

	int64_t shell_num_test;
	int64_t prim_num_test;
	int64_t ao_num_test;
	int64_t *nucleus_index_test;
	int64_t *nucleus_shell_num_test;
	int32_t *shell_ang_mom_test;
	int64_t *shell_prim_num_test;
	int64_t *shell_prim_index_test;
	double *shell_factor_test;
	double *exponent_test;
	double *coefficient_test;
	double *prim_factor_test;
	double *ao_factor_test;
	char typ_test;

	rc = qmckl_get_ao_basis_type_device(context, &typ_test);
	if (rc != QMCKL_SUCCESS)
		return 1;
	if (typ != typ_test)
		return 1;

	rc = qmckl_get_ao_basis_shell_num_device(context, &shell_num_test);
	if (rc != QMCKL_SUCCESS)
		return 1;
	if (shell_num != shell_num_test)
		return 1;

	rc = qmckl_get_ao_basis_prim_num_device(context, &prim_num_test);
	if (rc != QMCKL_SUCCESS)
		return 1;
	if (prim_num != prim_num_test)
		return 1;

	nucleus_index_test =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	rc = qmckl_get_ao_basis_nucleus_index_device(context, nucleus_index_test,
												 nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	bool wrong_val = false;

#pragma omp target parallel for is_device_ptr(nucleus_index_d,                 \
											  nucleus_index_test)
	for (int64_t i = 0; i < nucl_num; ++i) {
		if (nucleus_index_test[i] != nucleus_index_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, nucleus_index_test);
	if (wrong_val)
		return 1;

	nucleus_shell_num_test =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));
	rc = qmckl_get_ao_basis_nucleus_shell_num_device(
		context, nucleus_shell_num_test, nucl_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(nucleus_shell_num_d,             \
											  nucleus_shell_num_test)
	for (int64_t i = 0; i < nucl_num; ++i) {
		if (nucleus_shell_num_test[i] != nucleus_shell_num_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, nucleus_shell_num_test);
	if (wrong_val)
		return 1;

	shell_ang_mom_test =
		(int32_t *)qmckl_malloc_device(context, shell_num * sizeof(int32_t));
	rc = qmckl_get_ao_basis_shell_ang_mom_device(context, shell_ang_mom_test,
												 shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(shell_ang_mom_d,                 \
											  shell_ang_mom_test)
	for (int64_t i = 0; i < shell_num; ++i) {
		if (shell_ang_mom_test[i] != shell_ang_mom_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, shell_ang_mom_test);
	if (wrong_val)
		return 1;

	shell_factor_test =
		(double *)qmckl_malloc_device(context, shell_num * sizeof(double));
	rc = qmckl_get_ao_basis_shell_factor_device(context, shell_factor_test,
												shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(shell_factor_d, shell_factor_test)
	for (int64_t i = 0; i < shell_num; ++i) {
		if (shell_factor_test[i] != shell_factor_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, shell_factor_test);
	if (wrong_val)
		return 1;

	shell_prim_num_test =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	rc = qmckl_get_ao_basis_shell_prim_num_device(context, shell_prim_num_test,
												  shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

	shell_prim_index_test =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));
	rc = qmckl_get_ao_basis_shell_prim_index_device(
		context, shell_prim_index_test, shell_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(shell_prim_index_d,              \
											  shell_prim_index_test)
	for (int64_t i = 0; i < shell_num; ++i) {
		if (shell_prim_index_test[i] != shell_prim_index_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, shell_prim_index_test);
	if (wrong_val)
		return 1;

	exponent_test =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));
	rc = qmckl_get_ao_basis_exponent_device(context, exponent_test, prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(exponent_d, exponent_test)
	for (int64_t i = 0; i < prim_num; ++i) {
		if (exponent_test[i] != exponent_d[i])
			wrong_val = true;
		;
	}
	qmckl_free_device(context, exponent_test);
	if (wrong_val)
		return 1;

	coefficient_test =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));
	rc = qmckl_get_ao_basis_coefficient_device(context, coefficient_test,
											   prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(coefficient_d, coefficient_test)
	for (int64_t i = 0; i < prim_num; ++i) {
		if (coefficient_test[i] != coefficient_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, coefficient_test);
	if (wrong_val)
		return 1;

	prim_factor_test =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));
	rc = qmckl_get_ao_basis_prim_factor_device(context, prim_factor_test,
											   prim_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(prim_factor_d, prim_factor_test)
	for (int64_t i = 0; i < prim_num; ++i) {
		if (prim_factor_test[i] != prim_factor_d[i])
			wrong_val = true;
	}
	qmckl_free_device(context, prim_factor_test);
	if (wrong_val)
		return 1;

	rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num_test);
	if (ao_num != ao_num_test)
		return 1;

	ao_factor_test =
		(double *)qmckl_malloc_device(context, ao_num * sizeof(double));
	rc = qmckl_get_ao_basis_ao_factor_device(context, ao_factor_test, ao_num);
	if (rc != QMCKL_SUCCESS)
		return 1;

#pragma omp target parallel for is_device_ptr(ao_factor_d, ao_factor_test)
	for (int64_t i = 0; i < ao_num; ++i) {
		if (ao_factor_test[i] != ao_factor_d[i])
			wrong_val = true;
		;
	}
	qmckl_free_device(context, ao_factor_test);
	if (wrong_val)
		return 1;

		// Test ao_vgl  values
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num
#define elec_num chbrclf_elec_num

	double *elec_coord = &(chbrclf_elec_coord[0][0][0]);
	const int64_t point_num = elec_num;
	const double *elec_coord_d =
		qmckl_malloc_device(context, 3 * point_num * sizeof(double));

	qmckl_memcpy_H2D(context, elec_coord_d, elec_coord,
					 3 * point_num * sizeof(double));

	// TODO Fix this
	// if (!qmckl_electron_provided(context))
	//	return 1;

	rc = qmckl_set_point_device(context, 'N', point_num, elec_coord_d,
								point_num * 3);
	if (rc != QMCKL_SUCCESS)
		return 1;

	double *ao_vgl_d =
		qmckl_malloc_device(context, point_num * 5 * ao_num * sizeof(double));
	double *ao_vgl =
		qmckl_malloc_host(context, point_num * 5 * ao_num * sizeof(double));

	rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_d,
										  (int64_t)5 * point_num * ao_num);

	qmckl_memcpy_D2H(context, ao_vgl, ao_vgl_d,
					 point_num * 5 * ao_num * sizeof(double));
	if (rc != QMCKL_SUCCESS)
		return 1;

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
	printf(" ao_vgl ao_vgl[26][0][220] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 220)]);
	printf(" ao_vgl ao_vgl[26][1][220] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 220)]);
	printf(" ao_vgl ao_vgl[26][2][220] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 220)]);
	printf(" ao_vgl ao_vgl[26][3][220] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 220)]);
	printf(" ao_vgl ao_vgl[26][4][220] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 220)]);
	printf(" ao_vgl ao_vgl[26][0][221] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 221)]);
	printf(" ao_vgl ao_vgl[26][1][221] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 221)]);
	printf(" ao_vgl ao_vgl[26][2][221] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 221)]);
	printf(" ao_vgl ao_vgl[26][3][221] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 221)]);
	printf(" ao_vgl ao_vgl[26][4][221] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 221)]);
	printf(" ao_vgl ao_vgl[26][0][222] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 222)]);
	printf(" ao_vgl ao_vgl[26][1][222] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 222)]);
	printf(" ao_vgl ao_vgl[26][2][222] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 222)]);
	printf(" ao_vgl ao_vgl[26][3][222] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 222)]);
	printf(" ao_vgl ao_vgl[26][4][222] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 222)]);
	printf(" ao_vgl ao_vgl[26][0][223] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 223)]);
	printf(" ao_vgl ao_vgl[26][1][223] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 223)]);
	printf(" ao_vgl ao_vgl[26][2][223] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 223)]);
	printf(" ao_vgl ao_vgl[26][3][223] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 223)]);
	printf(" ao_vgl ao_vgl[26][4][223] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 223)]);
	printf(" ao_vgl ao_vgl[26][0][224] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 0, 224)]);
	printf(" ao_vgl ao_vgl[26][1][224] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 1, 224)]);
	printf(" ao_vgl ao_vgl[26][2][224] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 2, 224)]);
	printf(" ao_vgl ao_vgl[26][3][224] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 3, 224)]);
	printf(" ao_vgl ao_vgl[26][4][224] %25.15e\n",
		   ao_vgl[AO_VGL_ID(26, 4, 224)]);
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
	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 220)] - (1.516643537739178e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 220)] - (-7.725221462603871e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 220)] - (-6.507140835104833e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 220)] - (2.154644255710413e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 220)] - (6.365449359656352e-07)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 221)] - (-4.686370882518819e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 221)] - (2.387064067626827e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 221)] - (2.154644255710412e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 221)] - (-1.998731863512374e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 221)] - (-1.966899656441993e-07)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 222)] - (7.514816980753531e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 222)] - (-4.025889138635182e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 222)] - (-2.993372555126361e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 222)] - (1.067604670272904e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 222)] - (3.168199650002648e-07)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 223)] - (-4.021908374204471e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 223)] - (2.154644255710413e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 223)] - (1.725594944732276e-08)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 223)] - (-1.715339357718333e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 223)] - (-1.688020516893476e-07)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 0, 224)] - (7.175045873560788e-10)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 1, 224)] - (-3.843864637762753e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 2, 224)] - (-3.298857850451910e-09)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 3, 224)] - (-4.073047518790881e-10)) > 1.e-14)
		return 1;
	if (fabs(ao_vgl[AO_VGL_ID(26, 4, 224)] - (3.153244195820293e-08)) > 1.e-14)
		return 1;

	// TODO Fix this
	// rc = qmckl_context_destroy_device(context);
	// if (rc != QMCKL_SUCCESS)
	// 	return 1;

	return 0;
}
