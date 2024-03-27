#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "n2.h"
#include "../include/qmckl_gpu.h"
#include <omp.h>

int main() {

	qmckl_context_device context;

	if (omp_get_num_devices() <= 0) {
		printf("Error : no device found. Aborting execution\n");
		exit(1);
	}
	context = qmckl_context_create_device(omp_get_default_device());

	/* Reference input data */

	int64_t walk_num = n2_walk_num;
	int64_t elec_num = n2_elec_num;
	int64_t elec_up_num = n2_elec_up_num;
	int64_t elec_dn_num = n2_elec_dn_num;
	int64_t nucl_num = n2_nucl_num;
	double rescale_factor_ee = 0.6;

	double *rescale_factor_en =
		qmckl_malloc_device(context, 2 * sizeof(double));
#pragma omp target is_device_ptr(rescale_factor_en)
	{
		rescale_factor_en[0] = 0.6;
		rescale_factor_en[1] = 0.6;
	}

	// double *elec_coord = &(n2_elec_coord[0][0][0]);
	double *elec_coord = qmckl_malloc_device(context, sizeof(n2_elec_coord));
	qmckl_memcpy_H2D(context, elec_coord, n2_elec_coord, sizeof(n2_elec_coord));

	// double *nucl_charge = n2_charge;
	double *nucl_charge = qmckl_malloc_device(context, sizeof(n2_charge));
	qmckl_memcpy_H2D(context, nucl_charge, n2_charge, sizeof(n2_charge));

	// double *nucl_coord = &(n2_nucl_coord[0][0]);
	double *nucl_coord = qmckl_malloc_device(context, sizeof(n2_nucl_coord));
	qmckl_memcpy_H2D(context, nucl_coord, n2_nucl_coord, sizeof(n2_nucl_coord));

	int64_t size_max;

	/* Provide Electron data */

	qmckl_exit_code_device rc;

	rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);

	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord,
										 walk_num * 3 * elec_num);

	double *elec_coord2 =
		qmckl_malloc_device(context, walk_num * 3 * elec_num * sizeof(double));

	rc = qmckl_get_electron_coord_device(context, 'N', elec_coord2,
										 walk_num * 3 * elec_num);

	bool *wrongval_d = qmckl_malloc_device(context, sizeof(bool));
	bool *wrongval_h = malloc(sizeof(bool));

#pragma omp target is_device_ptr(wrongval_d)
	{ wrongval_d[0] = false; }
	wrongval_h[0] = false;

#pragma omp target is_device_ptr(elec_coord, elec_coord2, wrongval_d)
	{
		for (int64_t i = 0; i < 3 * elec_num; ++i) {
			if (elec_coord[i] != elec_coord2[i]) {
				wrongval_d[0] = true;
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	/* Provide Nucleus data */

	rc = qmckl_set_nucleus_num_device(context, nucl_num);

	double *nucl_coord2 =
		qmckl_malloc_device(context, 3 * nucl_num * sizeof(double));

	rc =
		qmckl_get_nucleus_coord_device(context, 'T', nucl_coord2, 3 * nucl_num);

	rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord, 3 * nucl_num);

	rc =
		qmckl_get_nucleus_coord_device(context, 'N', nucl_coord2, nucl_num * 3);

#pragma omp target is_device_ptr(nucl_coord, nucl_coord2, wrongval_d)
	{
		for (int64_t k = 0; k < 3; ++k) {
			for (int64_t i = 0; i < nucl_num; ++i) {
				if (nucl_coord[nucl_num * k + i] != nucl_coord2[3 * i + k]) {
					wrongval_d[0] = true;
				}
				if (wrongval_d[0]) {
				}
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	rc =
		qmckl_get_nucleus_coord_device(context, 'T', nucl_coord2, nucl_num * 3);
#pragma omp target is_device_ptr(nucl_coord, nucl_coord2, wrongval_d)
	{
		for (int64_t i = 0; i < 3 * nucl_num; ++i) {
			if (nucl_coord[i] != nucl_coord2[i]) {
				wrongval_d[0] = true;
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *nucl_charge2 = qmckl_malloc_device(context, sizeof(n2_charge));

	rc = qmckl_get_nucleus_charge_device(context, nucl_charge2, nucl_num);

	rc = qmckl_set_nucleus_charge_device(context, nucl_charge, nucl_num);

	rc = qmckl_get_nucleus_charge_device(context, nucl_charge2, nucl_num);
	for (int64_t i = 0; i < nucl_num; ++i) {
#pragma omp target is_device_ptr(nucl_charge, nucl_charge2, wrongval_d)
		{
			if (nucl_charge[i] != nucl_charge2[i]) {
				wrongval_d[0] = true;
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	int64_t type_nucl_num = n2_type_nucl_num;

	// int64_t *type_nucl_vector = &(n2_type_nucl_vector[0]);
	int64_t *type_nucl_vector =
		qmckl_malloc_device(context, sizeof(n2_type_nucl_vector));
	qmckl_memcpy_H2D(context, type_nucl_vector, n2_type_nucl_vector,
					 sizeof(n2_type_nucl_vector));

	int64_t aord_num = n2_aord_num;
	int64_t bord_num = n2_bord_num;
	int64_t cord_num = n2_cord_num;

	// double *a_vector = &(n2_aord_vector[0][0]);
	double *a_vector = qmckl_malloc_device(context, sizeof(n2_a_vector));
	qmckl_memcpy_H2D(context, a_vector, n2_a_vector, sizeof(n2_a_vector));

	// double *b_vector = &(n2_bord_vector[0]);
	double *b_vector = qmckl_malloc_device(context, sizeof(n2_b_vector));
	qmckl_memcpy_H2D(context, b_vector, n2_b_vector, sizeof(n2_b_vector));

	// double *c_vector = &(n2_cord_vector[0][0]);
	double *c_vector = qmckl_malloc_device(context, sizeof(n2_c_vector));
	qmckl_memcpy_H2D(context, c_vector, n2_c_vector, sizeof(n2_c_vector));

	int64_t dim_c_vector = 0;

	/* Set the data */
	rc = qmckl_set_jastrow_aord_num_device(context, aord_num);
	rc = qmckl_set_jastrow_bord_num_device(context, bord_num);
	rc = qmckl_set_jastrow_cord_num_device(context, cord_num);
	rc = qmckl_set_jastrow_type_nucl_num_device(context, type_nucl_num);
	rc = qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector,
												   nucl_num);
	rc = qmckl_set_jastrow_a_vector_device(context, a_vector,
										   (aord_num + 1) * type_nucl_num);
	rc = qmckl_set_jastrow_b_vector_device(context, b_vector, (bord_num + 1));
	rc = qmckl_get_jastrow_dim_c_vector_device(context, &dim_c_vector);
	rc = qmckl_set_jastrow_c_vector_device(context, c_vector,
										   dim_c_vector * type_nucl_num);

	double k_ee = 0.;

	double *k_en = qmckl_malloc_device(context, 2 * sizeof(double));
#pragma omp target is_device_ptr(k_en)
	{
		k_en[0] = 0.;
		k_en[1] = 0.;
	}

	rc = qmckl_set_jastrow_rescale_factor_en_device(context, rescale_factor_en,
													type_nucl_num);

	rc = qmckl_set_jastrow_rescale_factor_ee_device(context, rescale_factor_ee);

	rc = qmckl_get_jastrow_rescale_factor_ee_device(context, &k_ee);
	if (k_ee != rescale_factor_ee) {
		return 1;
	}

	rc = qmckl_get_jastrow_rescale_factor_en_device(context, k_en,
													type_nucl_num);

#pragma omp target is_device_ptr(k_en, rescale_factor_en, wrongval_d)
	{
		for (int i = 0; i < type_nucl_num; ++i) {
			if (fabs(k_en[i] - rescale_factor_en[i]) > 1e-12) {
				wrongval_d[0] = true;
			}
		}
	}
	if (wrongval_h[0]) {
		return 1;
	}

	double *asymp_jasb = qmckl_malloc_device(context, 2 * sizeof(double));

	// calculate asymp_jasb
	rc = qmckl_get_jastrow_asymp_jasb_device(context, asymp_jasb, 2);
#pragma omp target is_device_ptr(asymp_jasb, wrongval_d)
	{
		if (fabs(asymp_jasb[0] - 0.7115733522582638) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(asymp_jasb[1] - 1.043287918508297) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *factor_ee = qmckl_malloc_device(context, walk_num * sizeof(double));
	rc = qmckl_get_jastrow_factor_ee_device(context, factor_ee, walk_num);

	// calculate factor_ee
	rc = qmckl_get_jastrow_factor_ee_device(context, factor_ee, walk_num);
#pragma omp target is_device_ptr(factor_ee, wrongval_d)
	{
		if (fabs(factor_ee[0] + 16.83886184243964) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	// calculate factor_ee_deriv_e
	double *factor_ee_deriv_e =
		qmckl_malloc_device(context, walk_num * 4 * elec_num * sizeof(double));

	// check factor_ee_deriv_e
	rc = qmckl_get_jastrow_factor_ee_deriv_e_device(context, factor_ee_deriv_e,
													walk_num * 4 * elec_num);

#pragma omp target is_device_ptr(factor_ee_deriv_e, wrongval_d)
	{

		if (fabs(factor_ee_deriv_e[0 + 0 * elec_num + 0] +
				 0.39319353942687446) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_ee_deriv_e[0 + 1 * elec_num + 0] - 1.0535615450668214) >
			1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_ee_deriv_e[0 + 2 * elec_num + 0] +
				 0.39098406960784515) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_ee_deriv_e[0 + 3 * elec_num + 0] - 2.8650469630854483) >
			1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *ee_distance_rescaled = qmckl_malloc_device(
		context, walk_num * elec_num * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_ee_distance_rescaled_device(context,
													   ee_distance_rescaled);

#pragma omp target is_device_ptr(ee_distance_rescaled, wrongval_d)
	{
		// (e1,e2,w)
		// (0,0,0) == 0.
		if (ee_distance_rescaled[0] != 0.) {
			wrongval_d[0] = true;
		}

		// (1,0,0) == (0,1,0)
		if (ee_distance_rescaled[1] != ee_distance_rescaled[elec_num]) {
			wrongval_d[0] = true;
		}

		// value of (1,0,0)
		if (fabs(ee_distance_rescaled[1] - 0.6347507420688708) > 1.e-12) {
			wrongval_d[0] = true;
		}

		// (0,0,1) == 0.
		if (ee_distance_rescaled[5 * elec_num + 5] != 0.) {
			wrongval_d[0] = true;
		}

		// (1,0,1) == (0,1,1)
		if (ee_distance_rescaled[5 * elec_num + 6] !=
			ee_distance_rescaled[6 * elec_num + 5]) {
			wrongval_d[0] = true;
		}

		// value of (1,0,1)
		if (fabs(ee_distance_rescaled[5 * elec_num + 6] - 0.3941735387855409) >
			1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *ee_distance_rescaled_deriv_e = qmckl_malloc_device(
		context, 4 * walk_num * elec_num * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device(
		context, ee_distance_rescaled_deriv_e);

	double *factor_en = qmckl_malloc_device(context, walk_num * sizeof(double));
	rc = qmckl_get_jastrow_factor_en_device(context, factor_en, walk_num);

	// calculate factor_en

#pragma omp target is_device_ptr(factor_en, wrongval_d)
	{
		// BUG fails at 1.e-4
		if (fabs(22.781375792083587 - factor_en[0]) > 1.e-4) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	// calculate factor_en_deriv_e
	double *factor_en_deriv_e =
		qmckl_malloc_device(context, walk_num * 4 * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_factor_en_deriv_e_device(context, factor_en_deriv_e,
													walk_num * 4 * elec_num);

// check factor_en_deriv_e
#pragma omp target is_device_ptr(factor_en_deriv_e, wrongval_d)
	{
		if (fabs(factor_en_deriv_e[0 + 0 * elec_num + 0] -
				 0.19656663796630847) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_en_deriv_e[0 + 1 * elec_num + 0] + 0.3945140890522283) >
			1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_en_deriv_e[0 + 2 * elec_num + 0] - 0.5082964671286118) >
			1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(factor_en_deriv_e[0 + 3 * elec_num + 0] + 1.8409460670666289) >
			1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *en_distance_rescaled = qmckl_malloc_device(
		context, walk_num * nucl_num * elec_num * sizeof(double));

	rc = qmckl_get_electron_en_distance_rescaled_device(context,
														en_distance_rescaled);

	// BUG fails at 1.e-4
#pragma omp target is_device_ptr(en_distance_rescaled, wrongval_d)
	{
		// (e,n,w) in Fortran notation
		// (1,1,1)
		if (fabs(en_distance_rescaled[0 + 0 + 0] - 0.4942158656729477) >
			1.e-12) {
			wrongval_d[0] = true;
		}

		// (1,2,1)
		if (fabs(en_distance_rescaled[0 + 1 * elec_num + 0] -
				 1.2464137498005765) > 1.e-4) {
			wrongval_d[0] = true;
		}

		// (2,1,1)
		if (fabs(en_distance_rescaled[0 + 0 + 1] - 0.5248654474756858) >
			1.e-12) {
			wrongval_d[0] = true;
		}

		// (1,1,2)
		if (fabs(en_distance_rescaled[0 + 0 + 5] - 0.19529459944794733) >
			1.e-12) {
			wrongval_d[0] = true;
		}

		// (1,2,2)
		if (fabs(en_distance_rescaled[0 + 1 * elec_num + 5] -
				 1.2091967687767369) > 1.e-12) {
			wrongval_d[0] = true;
		}

		// (2,1,2)
		if (fabs(en_distance_rescaled[0 + 0 + 6] - 0.4726452953409436) >
			1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *en_distance_rescaled_deriv_e = qmckl_malloc_device(
		context, walk_num * 4 * nucl_num * elec_num * sizeof(double));

	rc = qmckl_get_electron_en_distance_rescaled_deriv_e_device(
		context, en_distance_rescaled_deriv_e);

	double *een_rescaled_e =
		qmckl_malloc_device(context, walk_num * (cord_num + 1) * elec_num *
										 elec_num * sizeof(double));
	rc = qmckl_get_jastrow_een_rescaled_e_device(context, een_rescaled_e,
												 elec_num * elec_num *
													 (cord_num + 1) * walk_num);

#pragma omp target is_device_ptr(een_rescaled_e, wrongval_d)
	{

		// value of (0,1,0,2)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								1 * elec_num * elec_num + 0 * elec_num + 2] -
				 0.2211015082992776) > 1.e-12) {
			wrongval_d[0] = true;
		}
		// value of (0,1,0,3)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								1 * elec_num * elec_num + 0 * elec_num + 3] -
				 0.2611178387068169) > 1.e-12) {
			wrongval_d[0] = true;
		}
		// value of (0,1,0,4)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								1 * elec_num * elec_num + 0 * elec_num + 4] -
				 0.0884012350763747) > 1.e-12) {
			wrongval_d[0] = true;
		}
		// value of (0,2,1,3)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								2 * elec_num * elec_num + 1 * elec_num + 3] -
				 0.1016685507354656) > 1.e-12) {
			wrongval_d[0] = true;
		}
		// value of (0,2,1,4)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								2 * elec_num * elec_num + 1 * elec_num + 4] -
				 0.0113118073246869) > 1.e-12) {
			wrongval_d[0] = true;
		}
		// value of (0,2,1,5)
		if (fabs(een_rescaled_e[0 * elec_num * elec_num * (cord_num + 1) +
								2 * elec_num * elec_num + 1 * elec_num + 5] -
				 0.5257156022077619) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *een_rescaled_e_deriv_e =
		qmckl_malloc_device(context, walk_num * (cord_num + 1) * elec_num * 4 *
										 sizeof(double) * elec_num);
	size_max = walk_num * (cord_num + 1) * elec_num * 4 * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_e_deriv_e_device(
		context, een_rescaled_e_deriv_e, size_max);

#pragma omp target is_device_ptr(een_rescaled_e_deriv_e, wrongval_d)
	{
		// value of (0,0,0,2,1)
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										1 * elec_num * 4 * elec_num +
										0 * elec_num * 4 + 0 * elec_num + 2] +
				 0.09831391870751387) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										1 * elec_num * 4 * elec_num +
										0 * elec_num * 4 + 0 * elec_num + 3] +
				 0.017204157459682526) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										1 * elec_num * 4 * elec_num +
										0 * elec_num * 4 + 0 * elec_num + 4] +
				 0.013345768421098641) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										2 * elec_num * 4 * elec_num +
										1 * elec_num * 4 + 0 * elec_num + 3] +
				 0.03733086358273962) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										2 * elec_num * 4 * elec_num +
										1 * elec_num * 4 + 0 * elec_num + 4] +
				 0.004922634822943517) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_e_deriv_e[0 * elec_num * 4 * elec_num *
											(cord_num + 1) +
										2 * elec_num * 4 * elec_num +
										1 * elec_num * 4 + 0 * elec_num + 5] +
				 0.5416751547830984) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *een_rescaled_n =
		qmckl_malloc_device(context, walk_num * (cord_num + 1) * nucl_num *
										 elec_num * sizeof(double));
	size_max = walk_num * (cord_num + 1) * nucl_num * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_n_device(context, een_rescaled_n,
												 size_max);

#pragma omp target is_device_ptr(een_rescaled_n, wrongval_d)
	{
		// value of (0,2,1)
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								1 * elec_num * nucl_num + 0 * elec_num + 2] -
				 0.2603169838750542) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								1 * elec_num * nucl_num + 0 * elec_num + 3] -
				 0.3016180139679065) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								1 * elec_num * nucl_num + 0 * elec_num + 4] -
				 0.10506023826192266) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								2 * elec_num * nucl_num + 1 * elec_num + 3] -
				 0.9267719759374164) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								2 * elec_num * nucl_num + 1 * elec_num + 4] -
				 0.11497585238132658) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n[0 * elec_num * nucl_num * (cord_num + 1) +
								2 * elec_num * nucl_num + 1 * elec_num + 5] -
				 0.07534033469115217) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *een_rescaled_n_deriv_e =
		qmckl_malloc_device(context, walk_num * (cord_num + 1) * nucl_num * 4 *
										 elec_num * sizeof(double));
	size_max = walk_num * (cord_num + 1) * nucl_num * 4 * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_n_deriv_e_device(
		context, een_rescaled_n_deriv_e, size_max);

// value of (0,2,1)
#pragma omp target is_device_ptr(een_rescaled_n_deriv_e, wrongval_d)
	{
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										1 * elec_num * 4 * nucl_num +
										0 * elec_num * 4 + 0 * elec_num + 2] +
				 0.11234061209936878) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										1 * elec_num * 4 * nucl_num +
										0 * elec_num * 4 + 0 * elec_num + 3] -
				 0.0004440109367151707) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										1 * elec_num * 4 * nucl_num +
										0 * elec_num * 4 + 0 * elec_num + 4] +
				 0.012868642597346566) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										2 * elec_num * 4 * nucl_num +
										1 * elec_num * 4 + 0 * elec_num + 3] -
				 0.08601122289922644) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										2 * elec_num * 4 * nucl_num +
										1 * elec_num * 4 + 0 * elec_num + 4] +
				 0.058681563677207206) > 1.e-12) {
			wrongval_d[0] = true;
		}
		if (fabs(een_rescaled_n_deriv_e[0 * elec_num * 4 * nucl_num *
											(cord_num + 1) +
										2 * elec_num * 4 * nucl_num +
										1 * elec_num * 4 + 0 * elec_num + 5] -
				 0.005359281880312882) > 1.e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *tmp_c =
		qmckl_malloc_device(context, walk_num * cord_num * (cord_num + 1) *
										 nucl_num * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_tmp_c_device(context, tmp_c);

	double *dtmp_c = qmckl_malloc_device(
		context, walk_num * cord_num * (cord_num + 1) * nucl_num * 4 *
					 elec_num * sizeof(double));
	rc = qmckl_get_jastrow_dtmp_c_device(context, dtmp_c);

#pragma omp target is_device_ptr(tmp_c, dtmp_c, wrongval_d)
	{
		if (fabs(tmp_c[0 + 0 + 1 * nucl_num * elec_num + 0 + 0] - 3.954384) >
			1e-6) {
			wrongval_d[0] = true;
		}

		if (fabs(dtmp_c[elec_num * 4 * nucl_num * (cord_num + 1)] -
				 3.278657e-01) > 1e-6) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *factor_een =
		qmckl_malloc_device(context, walk_num * sizeof(double));
	rc = qmckl_get_jastrow_factor_een_device(context, factor_een, walk_num);

#pragma omp target is_device_ptr(factor_een, wrongval_d)
	{
		if (fabs(factor_een[0] + 0.382580260174321) > 1e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	double *factor_een_deriv_e =
		qmckl_malloc_device(context, 4 * walk_num * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_factor_een_deriv_e_device(
		context, factor_een_deriv_e, 4 * walk_num * elec_num);

#pragma omp target is_device_ptr(factor_een_deriv_e, wrongval_d)
	{
		printf("%20.15e\n", factor_een_deriv_e[0 + 0 + 0]);
		if (fabs(factor_een_deriv_e[0 + 0 + 0] - 8.967809309100624e-02) >
			1e-12) {
			wrongval_d[0] = true;
		}

		printf("%20.15e\n",
			   factor_een_deriv_e[1 * elec_num * walk_num + 0 + 1]);
		if (fabs(factor_een_deriv_e[1 * elec_num * walk_num + 0 + 1] -
				 3.543090132452453e-02) > 1e-12) {
			wrongval_d[0] = true;
		}

		printf("%20.15e\n",
			   factor_een_deriv_e[2 * elec_num * walk_num + 0 + 2]);
		if (fabs(factor_een_deriv_e[2 * elec_num * walk_num + 0 + 2] -
				 8.996044894431991e-04) > 1e-12) {
			wrongval_d[0] = true;
		}

		printf("%20.15e\n",
			   factor_een_deriv_e[3 * elec_num * walk_num + 0 + 3]);
		if (fabs(factor_een_deriv_e[3 * elec_num * walk_num + 0 + 3] -
				 (-1.175028308456619e+00)) > 1e-12) {
			wrongval_d[0] = true;
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	printf("Total Jastrow value\n");

	rc = qmckl_get_jastrow_factor_ee_device(context, factor_ee, walk_num);
	rc = qmckl_get_jastrow_factor_en_device(context, factor_en, walk_num);
	rc = qmckl_get_jastrow_factor_een_device(context, factor_een, walk_num);
	double *total_j = qmckl_malloc_device(context, walk_num * sizeof(double));
	rc = qmckl_get_jastrow_value_device(context, total_j, walk_num);

#pragma omp target is_device_ptr(total_j, factor_ee, factor_en, factor_een,    \
								 wrongval_d)
	{
		for (int64_t i = 0; i < walk_num; ++i) {
			if (total_j[i] - exp(factor_ee[i] + factor_en[i] + factor_een[i]) >
				1.e-12) {
				wrongval_d[0] = true;
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	printf("Total Jastrow derivatives\n");

	rc = qmckl_get_jastrow_factor_ee_deriv_e_device(context, factor_ee_deriv_e,
													walk_num * elec_num * 4);

	rc = qmckl_get_jastrow_factor_en_deriv_e_device(context, factor_en_deriv_e,
													walk_num * elec_num * 4);

	rc = qmckl_get_jastrow_factor_een_deriv_e_device(
		context, factor_een_deriv_e, walk_num * elec_num * 4);

	double *total_j_deriv =
		qmckl_malloc_device(context, walk_num * 4 * elec_num * sizeof(double));
	rc = qmckl_get_jastrow_gl_device(context, total_j_deriv,
									 walk_num * elec_num * 4);

	rc = qmckl_get_jastrow_value_device(context, total_j, walk_num);

#pragma omp target is_device_ptr(total_j_deriv, total_j, factor_ee_deriv_e,    \
								 factor_en_deriv_e, factor_een_deriv_e,        \
								 wrongval_d)
	{
		for (int64_t k = 0; k < walk_num; ++k) {
			for (int64_t m = 0; m < 4; ++m) {
				for (int64_t e = 0; e < elec_num; ++e) {
					if (m < 3) { /* test only gradients */
						if (total_j_deriv[e + m * walk_num + k * elec_num * 4] /
									total_j[k] -
								(factor_ee_deriv_e[e + m * walk_num +
												   k * elec_num * 4] +
								 factor_en_deriv_e[e + m * walk_num +
												   k * elec_num * 4] +
								 factor_een_deriv_e[e + m * walk_num +
													k * elec_num * 4]) >
							1.e-12) {
							wrongval_d[0] = true;
						}
					}
				}
			}
		}
	}
	qmckl_memcpy_D2H(context, wrongval_h, wrongval_d, sizeof(bool));
	if (wrongval_h[0]) {
		return 1;
	}

	rc = qmckl_context_destroy_device(context);

	return 0;
}
