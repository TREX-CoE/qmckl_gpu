#pragma once

// This file contains prototypes for the Jastrow functions
// TODO Get up to date with the latest Jastrow rework

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_types.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_distance.h"

qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context_device context);
bool qmckl_jastrow_provided_device(qmckl_context_device context);

//**********
// SETTERS
//**********

qmckl_exit_code_device
qmckl_set_jastrow_rescale_factor_ee_device(qmckl_context_device context,
										   const double kappa_ee);
qmckl_exit_code_device
qmckl_set_jastrow_rescale_factor_en_device(qmckl_context_device context,
										   const double *kappa_en,
										   const int64_t size_max);
qmckl_exit_code_device
qmckl_set_jastrow_aord_num_device(qmckl_context_device context,
								  const int64_t aord_num);
qmckl_exit_code_device
qmckl_set_jastrow_bord_num_device(qmckl_context_device context,
								  const int64_t bord_num);
qmckl_exit_code_device
qmckl_set_jastrow_cord_num_device(qmckl_context_device context,
								  const int64_t cord_num);
qmckl_exit_code_device
qmckl_set_jastrow_type_nucl_num_device(qmckl_context_device context,
									   const int64_t type_nucl_num);
qmckl_exit_code_device
qmckl_set_jastrow_type_nucl_vector_device(qmckl_context_device context,
										  const int64_t *type_nucl_vector,
										  const int64_t nucl_num);
qmckl_exit_code_device
qmckl_set_jastrow_a_vector_device(qmckl_context_device context,
								  const double *a_vector,
								  const int64_t size_max);
qmckl_exit_code_device
qmckl_set_jastrow_b_vector_device(qmckl_context_device context,
								  const double *b_vector,
								  const int64_t size_max);
qmckl_exit_code_device
qmckl_set_jastrow_c_vector_device(qmckl_context_device context,
								  const double *c_vector,
								  const int64_t size_max);

//**********
// GETTERS (basic)
//**********

qmckl_exit_code_device
qmckl_get_jastrow_aord_num_device(qmckl_context_device context,
								  int64_t *const aord_num);
qmckl_exit_code_device
qmckl_get_jastrow_bord_num_device(qmckl_context_device context,
								  int64_t *const bord_num);
qmckl_exit_code_device
qmckl_get_jastrow_cord_num_device(qmckl_context_device context,
								  int64_t *const bord_num);
qmckl_exit_code_device
qmckl_get_jastrow_type_nucl_num_device(qmckl_context_device context,
									   int64_t *const type_nucl_num);
qmckl_exit_code_device
qmckl_get_jastrow_type_nucl_vector_device(qmckl_context_device context,
										  int64_t *const type_nucl_num,
										  const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_a_vector_device(qmckl_context_device context,
								  double *const a_vector,
								  const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_b_vector_device(qmckl_context_device context,
								  double *const b_vector,
								  const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_c_vector_device(qmckl_context_device context,
								  double *const c_vector,
								  const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_rescale_factor_ee_device(const qmckl_context_device context,
										   double *const rescale_factor_ee);
qmckl_exit_code_device
qmckl_get_jastrow_rescale_factor_en_device(const qmckl_context_device context,
										   double *const rescale_factor_en,
										   const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_dim_c_vector_device(qmckl_context_device context,
									  int64_t *const dim_c_vector);

qmckl_exit_code_device
qmckl_get_jastrow_asymp_jasb_device(qmckl_context_device context,
									double *const asymp_jasb,
									const int64_t size_max);

qmckl_exit_code_device
qmckl_get_jastrow_jasa_device(qmckl_context_device context,
							  double *const asymp_jasa, const int64_t size_max);

//**********
// COMPUTES
//**********

// Finalize

qmckl_exit_code_device qmckl_compute_jastrow_asymp_jasa_device(
	const qmckl_context_device context, const int64_t aord_num,
	const int64_t type_nucl_num, const double *a_vector,
	const double *rescale_factor_en, double *const asymp_jasa);

qmckl_exit_code_device qmckl_compute_jastrow_asymp_jasb_device(
	const qmckl_context_device context, const int64_t bord_num,
	const double *b_vector, const double rescale_factor_ee,
	double *const asymp_jasb);

qmckl_exit_code_device
qmckl_finalize_jastrow_device(qmckl_context_device context);

// Total Jastrow
qmckl_exit_code_device
qmckl_compute_jastrow_value_device(const qmckl_context_device context,
								   const int64_t walk_num, const double *f_ee,
								   const double *f_en, const double *f_een,
								   double *const value);
qmckl_exit_code_device qmckl_compute_jastrow_gl_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const double *value, const double *gl_ee,
	const double *gl_en, const double *gl_een, double *const gl);

// Electron/electron component
qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *asymp_jasb, double *const factor_ee);

qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *ee_distance_rescaled_deriv_e,
	double *const factor_ee_deriv_e);

// Electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_en_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *asymp_jasa, double *const factor_en);
qmckl_exit_code_device qmckl_compute_jastrow_factor_en_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *en_distance_rescaled_deriv_e,
	double *const factor_en_deriv_e);

// Electron/electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_een_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *tmp_c,
	const double *een_rescaled_n, double *const factor_een);
qmckl_exit_code_device qmckl_compute_jastrow_factor_een_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *tmp_c,
	const double *dtmp_c, const double *een_rescaled_n,
	const double *een_rescaled_n_deriv_e, double *const factor_een_deriv_e);

// Electron/electron/nucleus deriv
qmckl_exit_code_device
qmckl_compute_jastrow_factor_een_rescaled_e_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t cord_num,
	const double rescale_factor_ee, const double *coord_ee,
	const double *ee_distance, const double *een_rescaled_e,
	double *const een_rescaled_e_deriv_e);
qmckl_exit_code_device
qmckl_compute_jastrow_factor_een_rescaled_n_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const int64_t cord_num,
	const double *rescale_factor_en, const double *coord_ee,
	const double *coord_en, const double *en_distance,
	const double *een_rescaled_n, double *const een_rescaled_n_deriv_e);

// Distances
qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled);
qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_deriv_e_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled_deriv_e);
qmckl_exit_code_device qmckl_compute_en_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled);
qmckl_exit_code_device qmckl_compute_en_distance_rescaled_deriv_e_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled_deriv_e);

qmckl_exit_code_device qmckl_compute_een_rescaled_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t cord_num,
	const double rescale_factor_ee, const double *ee_distance,
	double *const een_rescaled_e);
qmckl_exit_code_device qmckl_compute_een_rescaled_n_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const int64_t cord_num,
	const double *rescale_factor_en, const double *en_distance,
	double *const een_rescaled_n);
qmckl_exit_code_device qmckl_compute_c_vector_full_device(
	const qmckl_context_device context, const int64_t nucl_num,
	const int64_t dim_c_vector, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const double *c_vector,
	double *const c_vector_full);
qmckl_exit_code_device qmckl_compute_lkpm_combined_index_device(
	const qmckl_context_device context, const int64_t cord_num,
	const int64_t dim_c_vector, int64_t *const lkpm_combined_index);
qmckl_exit_code_device
qmckl_compute_tmp_c_device(const qmckl_context_device context,
						   const int64_t cord_num, const int64_t elec_num,
						   const int64_t nucl_num, const int64_t walk_num,
						   const double *een_rescaled_e,
						   const double *een_rescaled_n, double *const tmp_c);
qmckl_exit_code_device
qmckl_compute_dtmp_c_device(const qmckl_context_device context,
							const int64_t cord_num, const int64_t elec_num,
							const int64_t nucl_num, const int64_t walk_num,
							const double *een_rescaled_e_deriv_e,
							const double *een_rescaled_n, double *const dtmp_c);

//**********
// PROVIDE
//**********

// Finalize
qmckl_exit_code_device
qmckl_provide_jastrow_asymp_jasa_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_jastrow_asymp_jasb_device(qmckl_context_device context);

// Total Jastrow
qmckl_exit_code_device
qmckl_provide_jastrow_value_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_gl_device(qmckl_context_device context);

// Electron/electron component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_ee_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_factor_ee_deriv_e_device(qmckl_context_device context);

// Electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_en_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_factor_en_deriv_e_device(qmckl_context_device context);

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_een_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_factor_een_deriv_e_device(qmckl_context_device context);

// Distances
qmckl_exit_code_device
qmckl_provide_ee_distance_rescaled_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_en_distance_rescaled_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_en_distance_rescaled_deriv_e_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_een_rescaled_e_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_een_rescaled_n_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_een_rescaled_n_deriv_e_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_een_rescaled_e_deriv_e_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_c_vector_full_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_lkpm_combined_index_device(qmckl_context_device context);
qmckl_exit_code_device qmckl_provide_tmp_c_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_dtmp_c_device(qmckl_context_device context);

// Misc
qmckl_exit_code_device
qmckl_compute_dim_c_vector_device(const qmckl_context_device context,
								  const int64_t cord_num,
								  int64_t *const dim_c_vector);

// TODO This is not defined in Jastrow in QMCkl : find the file and incorporate
// this in the GPU lib
qmckl_exit_code_device qmckl_provide_ee_distance_device(
	qmckl_context_device context); // Required by provide_een_rescaled_e
qmckl_exit_code_device qmckl_provide_en_distance_device(
	qmckl_context_device context); // Required by provide_een_rescaled_n

qmckl_exit_code_device
qmckl_provide_ee_distance_rescaled_deriv_e_device(qmckl_context_device context);

//**********
// GETTERS (for computes)
//**********

// Total Jastrow
qmckl_exit_code_device
qmckl_get_jastrow_value_device(qmckl_context_device context,
							   double *const value, const int64_t size_max);
qmckl_exit_code_device qmckl_get_jastrow_gl_device(qmckl_context_device context,
												   double *const gl,
												   const int64_t size_max);

// Electron/electron component
qmckl_exit_code_device
qmckl_get_jastrow_factor_ee_device(qmckl_context_device context,
								   double *const factor_ee,
								   const int64_t size_max);

// Electron/nucleus component
qmckl_exit_code_device
qmckl_get_jastrow_factor_en_device(qmckl_context_device context,
								   double *const factor_en,
								   const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_factor_en_deriv_e_device(qmckl_context_device context,
										   double *const factor_en_deriv_e,
										   const int64_t size_max);

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_get_jastrow_factor_een_device(qmckl_context_device context,
									double *const factor_een,
									const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_factor_een_deriv_e_device(qmckl_context_device context,
											double *const factor_een_deriv_e,
											const int64_t size_max);

// Distances
qmckl_exit_code_device
qmckl_get_jastrow_ee_distance_rescaled_device(qmckl_context_device context,
											  double *const distance_rescaled);

qmckl_exit_code_device qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device(
	qmckl_context_device context, double *const distance_rescaled_deriv_e);

qmckl_exit_code_device
qmckl_get_electron_en_distance_rescaled_device(qmckl_context_device context,
											   double *distance_rescaled);
qmckl_exit_code_device qmckl_get_electron_en_distance_rescaled_deriv_e_device(
	qmckl_context_device context, double *distance_rescaled_deriv_e);

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_e_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_e_deriv_e_device(qmckl_context_device context,
												double *const distance_rescaled,
												const int64_t size_max);

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_n_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_n_deriv_e_device(qmckl_context_device context,
												double *const distance_rescaled,
												const int64_t size_max);

qmckl_exit_code_device
qmckl_get_jastrow_tmp_c_device(qmckl_context_device context,
							   double *const tmp_c);
qmckl_exit_code_device
qmckl_get_jastrow_dtmp_c_device(qmckl_context_device context,
								double *const dtmp_c);

qmckl_exit_code_device
qmckl_get_jastrow_factor_ee_deriv_e_device(qmckl_context_device context,
										   double *const factor_ee_deriv_e,
										   const int64_t size_max);
