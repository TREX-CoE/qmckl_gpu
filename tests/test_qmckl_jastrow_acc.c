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
#include <openacc.h>

int main() {

	qmckl_context_device context;
	acc_device_t device_type = acc_get_device_type();
	if (acc_get_num_devices(device_type) <= 0) {
		printf("Error : no device found. Aborting execution\n");
		exit(1);
	}
	acc_init(device_type);
	context = qmckl_context_create_device(0);

	/* Reference input data */

	int64_t walk_num = n2_walk_num;
	int64_t elec_num = n2_elec_num;
	int64_t elec_up_num = n2_elec_up_num;
	int64_t elec_dn_num = n2_elec_dn_num;
	int64_t nucl_num = n2_nucl_num;
	double rescale_factor_ee = 1.0;
	double rescale_factor_en[2] = {1.0, 1.0};
	double *elec_coord = &(n2_elec_coord[0][0][0]);

	const double *nucl_charge = n2_charge;
	double *nucl_coord = &(n2_nucl_coord[0][0]);
	int64_t size_max;

	/* Provide Electron data */

	qmckl_exit_code_device rc;

	assert(!qmckl_electron_provided_device(context));

	rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	assert(qmckl_electron_provided_device(context));

	rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord,
										 walk_num * 3 * elec_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double elec_coord2[walk_num * 3 * elec_num];

	rc = context, qmckl_get_electron_coord_device(context, 'N', elec_coord2,
												  walk_num * 3 * elec_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	for (int64_t i = 0; i < 3 * elec_num; ++i) {
		assert(elec_coord[i] == elec_coord2[i]);
	}

	/* Provide Nucleus data */

	assert(!qmckl_nucleus_provided_device(context));

	rc = qmckl_set_nucleus_num_device(context, nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	assert(!qmckl_nucleus_provided_device(context));

	double nucl_coord2[3 * nucl_num];

	rc =
		qmckl_get_nucleus_coord_device(context, 'T', nucl_coord2, 3 * nucl_num);
	assert(rc == QMCKL_NOT_PROVIDED_DEVICE);

	rc = qmckl_set_nucleus_coord_device(context, 'T', &(nucl_coord[0]), 3 * nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	assert(!qmckl_nucleus_provided_device(context));

	rc =
		qmckl_get_nucleus_coord_device(context, 'N', nucl_coord2, nucl_num * 3);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	for (int64_t k = 0; k < 3; ++k) {
		for (int64_t i = 0; i < nucl_num; ++i) {
			assert(nucl_coord[nucl_num * k + i] == nucl_coord2[3 * i + k]);
		}
	}

	rc =
		qmckl_get_nucleus_coord_device(context, 'T', nucl_coord2, nucl_num * 3);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	for (int64_t i = 0; i < 3 * nucl_num; ++i) {
		assert(nucl_coord[i] == nucl_coord2[i]);
	}

	double nucl_charge2[nucl_num];

	rc = qmckl_get_nucleus_charge_device(context, nucl_charge2, nucl_num);
	assert(rc == QMCKL_NOT_PROVIDED_DEVICE);

	rc = qmckl_set_nucleus_charge_device(context, nucl_charge, nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_nucleus_charge_device(context, nucl_charge2, nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	for (int64_t i = 0; i < nucl_num; ++i) {
		assert(nucl_charge[i] == nucl_charge2[i]);
	}
	assert(qmckl_nucleus_provided_device(context));

	assert(qmckl_electron_provided_device(context));

	int64_t type_nucl_num = n2_type_nucl_num;
	int64_t *type_nucl_vector = &(n2_type_nucl_vector[0]);
	int64_t aord_num = n2_aord_num;
	int64_t bord_num = n2_bord_num;
	int64_t cord_num = n2_cord_num;
	double *a_vector = &(n2_aord_vector[0][0]);
	double *b_vector = &(n2_bord_vector[0]);
	double *c_vector = &(n2_cord_vector[0][0]);
	int64_t dim_c_vector = 0;

	assert(!qmckl_jastrow_provided_device(context));

	/* Set the data */
	rc = qmckl_set_jastrow_aord_num_device(context, aord_num);
	rc = qmckl_set_jastrow_bord_num_device(context, bord_num);
	rc = qmckl_set_jastrow_cord_num_device(context, cord_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_set_jastrow_type_nucl_num_device(context, type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector,
												  nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_set_jastrow_a_vector_device(context, a_vector,
										  (aord_num + 1) * type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_set_jastrow_b_vector_device(context, b_vector, (bord_num + 1));
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_get_jastrow_dim_c_vector_device(context, &dim_c_vector);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	rc = qmckl_set_jastrow_c_vector_device(context, c_vector,
										   dim_c_vector * type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double k_ee = 0.;
	double k_en[2] = {0., 0.};
	rc = qmckl_set_jastrow_rescale_factor_en_device(
		context, rescale_factor_en, type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_set_jastrow_rescale_factor_ee_device(context,
														  rescale_factor_ee);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_rescale_factor_ee_device(context, &k_ee);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	assert(k_ee == rescale_factor_ee);

	rc = qmckl_get_jastrow_rescale_factor_en_device(context, &(k_en[0]),
														  type_nucl_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);
	for (int i = 0; i < type_nucl_num; ++i) {
		assert(k_en[i] == rescale_factor_en[i]);
	}

	/* Check if Jastrow is properly initialized */

	assert(qmckl_jastrow_provided_device(context));

	double asymp_jasb[2];
	rc = qmckl_get_jastrow_asymp_jasb_device(context, asymp_jasb, 2);

	// calculate asymp_jasb
	assert(fabs(asymp_jasb[0] - 0.5323750557252571) < 1.e-12);
	assert(fabs(asymp_jasb[1] - 0.31567342786262853) < 1.e-12);

	/* Check if Jastrow is properly initialized */

	assert(qmckl_jastrow_provided_device(context));

	double factor_ee[walk_num];
	rc = qmckl_get_jastrow_factor_ee(context, factor_ee, walk_num);

	// calculate factor_ee
	printf("%e\n%e\n\n", factor_ee[0], -4.282760865958113);
	assert(fabs(factor_ee[0] + 4.282760865958113) < 1.e-12);

	/* Check if Jastrow is properly initialized */

	assert(qmckl_jastrow_provided_device(context));

	// calculate factor_ee_deriv_e
	double factor_ee_deriv_e[walk_num][4][elec_num];
	rc = qmckl_get_jastrow_factor_ee_deriv_e_device(
		context, &(factor_ee_deriv_e[0][0][0]), walk_num * 4 * elec_num);

	// check factor_ee_deriv_e
	assert(fabs(factor_ee_deriv_e[0][0][0] - 0.16364894652107934) < 1.e-12);
	assert(fabs(factor_ee_deriv_e[0][1][0] + 0.6927548119830084) < 1.e-12);
	assert(fabs(factor_ee_deriv_e[0][2][0] - 0.073267755223968) < 1.e-12);
	assert(fabs(factor_ee_deriv_e[0][3][0] - 1.5111672803213185) < 1.e-12);

	assert(qmckl_electron_provided_device(context));

	double ee_distance_rescaled[walk_num * elec_num * elec_num];
	rc = qmckl_get_jastrow_ee_distance_rescaled_device(context,
													   ee_distance_rescaled);

	// (e1,e2,w)
	// (0,0,0) == 0.
	assert(ee_distance_rescaled[0] == 0.);

	// (1,0,0) == (0,1,0)
	assert(ee_distance_rescaled[1] == ee_distance_rescaled[elec_num]);

	// value of (1,0,0)
	assert(fabs(ee_distance_rescaled[1] - 0.5502278003524018) < 1.e-12);

	// (0,0,1) == 0.
	assert(ee_distance_rescaled[5 * elec_num + 5] == 0.);

	// (1,0,1) == (0,1,1)
	assert(ee_distance_rescaled[5 * elec_num + 6] ==
		   ee_distance_rescaled[6 * elec_num + 5]);

	// value of (1,0,1)
	assert(fabs(ee_distance_rescaled[5 * elec_num + 6] - 0.3622098222364193) <
		   1.e-12);

	assert(qmckl_electron_provided_device(context));

	double ee_distance_rescaled_deriv_e[4 * walk_num * elec_num * elec_num];
	rc = qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device(
		context, ee_distance_rescaled_deriv_e);

	// TODO: Get exact values
	//// (e1,e2,w)
	//// (0,0,0) == 0.
	// assert(ee_distance[0] == 0.);
	//
	//// (1,0,0) == (0,1,0)
	// assert(ee_distance[1] == ee_distance[elec_num]);
	//
	//// value of (1,0,0)
	// assert(fabs(ee_distance[1]-7.152322512964209) < 1.e-12);
	//
	//// (0,0,1) == 0.
	// assert(ee_distance[elec_num*elec_num] == 0.);
	//
	//// (1,0,1) == (0,1,1)
	// assert(ee_distance[elec_num*elec_num+1] ==
	// ee_distance[elec_num*elec_num+elec_num]);
	//
	//// value of (1,0,1)
	// assert(fabs(ee_distance[elec_num*elec_num+1]-6.5517646321055665)
	// < 1.e-12);

	double asymp_jasa[2];
	rc =
		qmckl_get_jastrow_asymp_jasa_device(context, asymp_jasa, type_nucl_num);

	// calculate asymp_jasb
	printf("%e %e\n", asymp_jasa[0], -0.548554);
	assert(fabs(-0.548554 - asymp_jasa[0]) < 1.e-12);

	/* Check if Jastrow is properly initialized */
	assert(qmckl_jastrow_provided_device(context));

	double factor_en[walk_num];
	rc = qmckl_get_jastrow_factor_en_device(context, factor_en, walk_num);

	// calculate factor_en
	assert(fabs(5.1052574308112755 - factor_en[0]) < 1.e-12);

	/* Check if Jastrow is properly initialized */
	assert(qmckl_jastrow_provided_device(context));

	// calculate factor_en_deriv_e
	double factor_en_deriv_e[walk_num][4][elec_num];
	rc = qmckl_get_jastrow_factor_en_deriv_e(
		context, &(factor_en_deriv_e[0][0][0]), walk_num * 4 * elec_num);

	// check factor_en_deriv_e
	assert(fabs(factor_en_deriv_e[0][0][0] - 0.11609919541763383) < 1.e-12);
	assert(fabs(factor_en_deriv_e[0][1][0] + 0.23301394780804574) < 1.e-12);
	assert(fabs(factor_en_deriv_e[0][2][0] - 0.17548337641865783) < 1.e-12);
	assert(fabs(factor_en_deriv_e[0][3][0] + 0.9667363412285741) < 1.e-12);

	assert(qmckl_electron_provided_device(context));
	assert(qmckl_nucleus_provided_device(context));

	double en_distance_rescaled[walk_num][nucl_num][elec_num];

	rc = qmckl_get_electron_en_distance_rescaled_device(
		context, &(en_distance_rescaled[0][0][0]));
	assert(rc == QMCKL_SUCCESS_DEVICE);

	// (e,n,w) in Fortran notation
	// (1,1,1)
	assert(fabs(en_distance_rescaled[0][0][0] - 0.4435709484118112) < 1.e-12);

	// (1,2,1)
	assert(fabs(en_distance_rescaled[0][1][0] - 0.8993601506374442) < 1.e-12);

	// (2,1,1)
	assert(fabs(en_distance_rescaled[0][0][1] - 0.46760219699910477) < 1.e-12);

	// (1,1,2)
	assert(fabs(en_distance_rescaled[0][0][5] - 0.1875631834682101) < 1.e-12);

	// (1,2,2)
	assert(fabs(en_distance_rescaled[0][1][5] - 0.8840716589810682) < 1.e-12);

	// (2,1,2)
	assert(fabs(en_distance_rescaled[0][0][6] - 0.42640469987268914) < 1.e-12);

	assert(qmckl_electron_provided_device(context));

	assert(qmckl_nucleus_provided_device(context));

	double en_distance_rescaled_deriv_e[walk_num][4][nucl_num][elec_num];

	rc = qmckl_get_electron_en_distance_rescaled_deriv_e_device(
		context, &(en_distance_rescaled_deriv_e[0][0][0][0]));
	assert(rc == QMCKL_SUCCESS_DEVICE);

	// TODO: check exact values
	//// (e,n,w) in Fortran notation
	//// (1,1,1)
	// assert(fabs(en_distance_rescaled[0][0][0] - 7.546738741619978) < 1.e-12);
	//
	//// (1,2,1)
	// assert(fabs(en_distance_rescaled[0][1][0] - 8.77102435246984) < 1.e-12);
	//
	//// (2,1,1)
	// assert(fabs(en_distance_rescaled[0][0][1] - 3.698922010513608) < 1.e-12);
	//
	//// (1,1,2)
	// assert(fabs(en_distance_rescaled[1][0][0] - 5.824059436060509) < 1.e-12);
	//
	//// (1,2,2)
	// assert(fabs(en_distance_rescaled[1][1][0] - 7.080482110317645) < 1.e-12);
	//
	//// (2,1,2)
	// assert(fabs(en_distance_rescaled[1][0][1] - 3.1804527583077356)
	// < 1.e-12);

	assert(qmckl_electron_provided_device(context));

	double een_rescaled_e[walk_num][(cord_num + 1)][elec_num][elec_num];
	rc = qmckl_get_jastrow_een_rescaled_e_device(
		context, &(een_rescaled_e[0][0][0][0]),
		elec_num * elec_num * (cord_num + 1) * walk_num);

	// value of (0,2,1)
	assert(fabs(een_rescaled_e[0][1][0][2] - 0.08084493981483197) < 1.e-12);
	assert(fabs(een_rescaled_e[0][1][0][3] - 0.1066745707571846) < 1.e-12);
	assert(fabs(een_rescaled_e[0][1][0][4] - 0.01754273169464735) < 1.e-12);
	assert(fabs(een_rescaled_e[0][2][1][3] - 0.02214680362033448) < 1.e-12);
	assert(fabs(een_rescaled_e[0][2][1][4] - 0.0005700154999202759) < 1.e-12);
	assert(fabs(een_rescaled_e[0][2][1][5] - 0.3424402276009091) < 1.e-12);

	double een_rescaled_e_deriv_e[walk_num][(cord_num + 1)][elec_num][4]
								 [elec_num];
	size_max = walk_num * (cord_num + 1) * elec_num * 4 * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_e_deriv_e_device(
		context, &(een_rescaled_e_deriv_e[0][0][0][0][0]), size_max);

	// value of (0,0,0,2,1)
	assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][2] + 0.05991352796887283) <
		   1.e-12);
	assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][3] + 0.011714035071545248) <
		   1.e-12);
	assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][4] + 0.00441398875758468) <
		   1.e-12);
	assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][3] + 0.013553180060167595) <
		   1.e-12);
	assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][4] +
				0.00041342909359870457) < 1.e-12);
	assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][5] + 0.5880599146214673) <
		   1.e-12);

	assert(qmckl_electron_provided_device(context));

	double een_rescaled_n[walk_num][(cord_num + 1)][nucl_num][elec_num];
	size_max = walk_num * (cord_num + 1) * nucl_num * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_n_device(
		context, &(een_rescaled_n[0][0][0][0]), size_max);

	// value of (0,2,1)
	assert(fabs(een_rescaled_n[0][1][0][2] - 0.10612983920006765) < 1.e-12);
	assert(fabs(een_rescaled_n[0][1][0][3] - 0.135652809635553) < 1.e-12);
	assert(fabs(een_rescaled_n[0][1][0][4] - 0.023391817607642338) < 1.e-12);
	assert(fabs(een_rescaled_n[0][2][1][3] - 0.880957224822116) < 1.e-12);
	assert(fabs(een_rescaled_n[0][2][1][4] - 0.027185942659395074) < 1.e-12);
	assert(fabs(een_rescaled_n[0][2][1][5] - 0.01343938025140174) < 1.e-12);

	assert(qmckl_electron_provided_device(context));

	double een_rescaled_n_deriv_e[walk_num][(cord_num + 1)][nucl_num][4]
								 [elec_num];
	size_max = walk_num * (cord_num + 1) * nucl_num * 4 * elec_num;
	rc = qmckl_get_jastrow_een_rescaled_n_deriv_e_device(
		context, &(een_rescaled_n_deriv_e[0][0][0][0][0]), size_max);

	// value of (0,2,1)
	assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][2] + 0.07633444246999128) <
		   1.e-12);
	assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][3] -
				0.00033282346259738276) < 1.e-12);
	assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][4] + 0.004775370547333061) <
		   1.e-12);
	assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][3] - 0.1362654644223866) <
		   1.e-12);
	assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][4] + 0.0231253431662794) <
		   1.e-12);
	assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][5] - 0.001593334817691633) <
		   1.e-12);

	assert(qmckl_electron_provided_device(context));

	double tmp_c[walk_num][cord_num][cord_num + 1][nucl_num][elec_num];
	rc = qmckl_get_jastrow_tmp_c_device(context, &(tmp_c[0][0][0][0][0]));

	double dtmp_c[walk_num][cord_num][cord_num + 1][nucl_num][4][elec_num];
	rc = qmckl_get_jastrow_dtmp_c_device(context, &(dtmp_c[0][0][0][0][0][0]));

	printf("%e\n%e\n", tmp_c[0][0][1][0][0], 2.7083473948352403);
	assert(fabs(tmp_c[0][0][1][0][0] - 2.7083473948352403) < 1e-12);

	printf("%e\n%e\n", dtmp_c[0][1][0][0][0][0], 0.237440520852232);
	assert(fabs(dtmp_c[0][1][0][0][0][0] - 0.237440520852232) < 1e-12);

	/* Check if Jastrow is properly initialized */
	assert(qmckl_jastrow_provided_device(context));

	double factor_een[walk_num];
	rc = qmckl_get_jastrow_factor_een_device(context, &(factor_een[0]),
											 walk_num);

	assert(fabs(factor_een[0] + 0.37407972141304213) < 1e-12);

	/* Check if Jastrow is properly initialized */
	assert(qmckl_jastrow_provided_device(context));

	double factor_een_deriv_e[4][walk_num][elec_num];
	rc = qmckl_get_jastrow_factor_een_deriv_e_device(
		context, &(factor_een_deriv_e[0][0][0]), 4 * walk_num * elec_num);

	printf("%20.15e\n", factor_een_deriv_e[0][0][0]);
	assert(fabs(factor_een_deriv_e[0][0][0] - (-5.481671107220383e-04)) <
		   1e-12);

	printf("%20.15e\n", factor_een_deriv_e[1][0][1]);
	assert(fabs(factor_een_deriv_e[1][0][1] - (-5.402107832095666e-02)) <
		   1e-12);

	printf("%20.15e\n", factor_een_deriv_e[2][0][2]);
	assert(fabs(factor_een_deriv_e[2][0][2] - (-1.648945927082279e-01)) <
		   1e-12);

	printf("%20.15e\n", factor_een_deriv_e[3][0][3]);
	assert(fabs(factor_een_deriv_e[3][0][3] - (-1.269746119491287e+00)) <
		   1e-12);

	printf("Total Jastrow value\n");
	/* Check if Jastrow is properly initialized */
	assert(qmckl_jastrow_provided_device(context));

	rc = qmckl_get_jastrow_factor_ee_device(context, &(factor_ee[0]), walk_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_factor_en_device(context, &(factor_en[0]), walk_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_factor_een_device(context, &(factor_een[0]),
											 walk_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double total_j[walk_num];
	rc = qmckl_get_jastrow_value_device(context, &(total_j[0]), walk_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	for (int64_t i = 0; i < walk_num; ++i) {
		assert(total_j[i] - exp(factor_ee[i] + factor_en[i] + factor_een[i]) <
			   1.e-12);
	}

	printf("Total Jastrow derivatives\n");

	/* Check if Jastrow is properly initialized */

	assert(qmckl_jastrow_provided_device(context));

	rc = qmckl_get_factor_ee_deriv_e_device(
		context, &(factor_ee_deriv_e[0][0][0]), walk_num * elec_num * 4);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_factor_en_deriv_e_device(
		context, &(factor_en_deriv_e[0][0][0]), walk_num * elec_num * 4);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_factor_een_deriv_e_device(
		context, &(factor_een_deriv_e[0][0][0]), walk_num * elec_num * 4);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	double total_j_deriv[walk_num][4][elec_num];
	rc = qmckl_get_jastrow_gl_device(context, &(total_j_deriv[0][0][0]),
										   walk_num * elec_num * 4);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	rc = qmckl_get_jastrow_value_device(context, &(total_j[0]), walk_num);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	for (int64_t k = 0; k < walk_num; ++k) {
		for (int64_t m = 0; m < 4; ++m) {
			for (int64_t e = 0; e < elec_num; ++e) {
				if (m < 3) { /* test only gradients */
					assert(total_j_deriv[k][m][e] / total_j[k] -
							   (factor_ee_deriv_e[k][m][e] +
								factor_en_deriv_e[k][m][e] +
								factor_een_deriv_e[k][m][e]) <
						   1.e-12);
				}
			}
		}
	}

	rc = qmckl_context_destroy_device(context);
	assert(rc == QMCKL_SUCCESS_DEVICE);

	return 0;
}
