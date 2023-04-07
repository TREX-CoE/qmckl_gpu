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

/* Init func */

qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context_device context);

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

//**********
// COMPUTES
//**********

// Finalize
qmckl_exit_code_device
qmckl_finalize_jastrow_device(qmckl_context_device context);

// Total Jastrow
qmckl_exit_code_device
qmckl_compute_jastrow_value_device(const qmckl_context_device context,
								   const int64_t walk_num, const double *f_ee,
								   const double *f_en, const double *f_een,
								   double *const value);

// Electron/electron component
qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *asymp_jasb, double *const factor_ee);

// Electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_en_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *asymp_jasa, double *const factor_en);

// Electron/electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_een_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *een_rescaled_e,
	const double *een_rescaled_n, double *const factor_een);

// Distances
qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled);
qmckl_exit_code_device qmckl_compute_en_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled);

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

//**********
// PROVIDE
//**********

// Total Jastrow
qmckl_exit_code_device
qmckl_provide_jastrow_value_device(qmckl_context_device context);

// Electron/electron component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_ee_device(qmckl_context_device context);

// Electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_en_device(qmckl_context_device context);

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_provide_jastrow_factor_een_device(qmckl_context_device context);

// Distances
qmckl_exit_code_device
qmckl_provide_ee_distance_rescaled_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_en_distance_rescaled_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_een_rescaled_e_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_een_rescaled_n_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_jastrow_c_vector_full_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_provide_lkpm_combined_index_device(qmckl_context_device context);
qmckl_exit_code_device qmckl_provide_tmp_c_device(qmckl_context_device context);

// TODO This is not defined in Jastrow in QMCkl : find the file and incorporate
// this in the GPU lib
qmckl_exit_code_device qmckl_provide_ee_distance_device(
	qmckl_context_device context); // Required by provide_een_rescaled_e
qmckl_exit_code_device qmckl_provide_en_distance_device(
	qmckl_context_device context); // Required by provide_een_rescaled_n

//**********
// GETTERS (for computes)
//**********

// Total Jastrow
qmckl_exit_code_device
qmckl_get_jastrow_value_device(qmckl_context_device context,
							   double *const value, const int64_t size_max);

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

// Electron/electron/nucleus component
qmckl_exit_code_device
qmckl_get_jastrow_factor_een_device(qmckl_context_device context,
									double *const factor_een,
									const int64_t size_max);

// Distances
qmckl_exit_code_device
qmckl_get_jastrow_ee_distance_rescaled_device(qmckl_context_device context,
											  double *const distance_rescaled);
qmckl_exit_code_device
qmckl_get_electron_en_distance_rescaled_device(qmckl_context_device context,
											   double *distance_rescaled);

qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_e_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_een_rescaled_n_device(qmckl_context_device context,
										double *const distance_rescaled,
										const int64_t size_max);
qmckl_exit_code_device
qmckl_get_jastrow_tmp_c_device(qmckl_context_device context,
							   double *const tmp_c);
