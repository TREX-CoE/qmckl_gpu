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

/* Init func */

qmckl_exit_code qmckl_init_jastrow_device(qmckl_context context);

//**********
// SETTERS
//**********

qmckl_exit_code
qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_context context,
										  const double kappa_ee);
qmckl_exit_code qmckl_set_jastrow_champ_rescale_factor_en(
	qmckl_context context, const double *kappa_en, const int64_t size_max);
qmckl_exit_code qmckl_set_jastrow_champ_aord_num(qmckl_context context,
												 const int64_t aord_num);
qmckl_exit_code qmckl_set_jastrow_champ_bord_num(qmckl_context context,
												 const int64_t bord_num);
qmckl_exit_code qmckl_set_jastrow_champ_cord_num(qmckl_context context,
												 const int64_t cord_num);
qmckl_exit_code
qmckl_set_jastrow_champ_type_nucl_num(qmckl_context context,
									  const int64_t type_nucl_num);
qmckl_exit_code
qmckl_set_jastrow_champ_type_nucl_vector(qmckl_context context,
										 const int64_t *type_nucl_vector,
										 const int64_t nucl_num);
qmckl_exit_code qmckl_set_jastrow_champ_a_vector(qmckl_context context,
												 const double *a_vector,
												 const int64_t size_max);
qmckl_exit_code qmckl_set_jastrow_champ_b_vector(qmckl_context context,
												 const double *b_vector,
												 const int64_t size_max);
qmckl_exit_code qmckl_set_jastrow_champ_c_vector(qmckl_context context,
												 const double *c_vector,
												 const int64_t size_max);

//**********
// GETTERS (basic)
//**********

qmckl_exit_code qmckl_get_jastrow_champ_aord_num(qmckl_context context,
												 int64_t *const aord_num);
qmckl_exit_code qmckl_get_jastrow_champ_bord_num(qmckl_context context,
												 int64_t *const bord_num);
qmckl_exit_code qmckl_get_jastrow_champ_cord_num(qmckl_context context,
												 int64_t *const bord_num);
qmckl_exit_code
qmckl_get_jastrow_champ_type_nucl_num(qmckl_context context,
									  int64_t *const type_nucl_num);
qmckl_exit_code
qmckl_get_jastrow_champ_type_nucl_vector(qmckl_context context,
										 int64_t *const type_nucl_num,
										 const int64_t size_max);
qmckl_exit_code qmckl_get_jastrow_champ_a_vector(qmckl_context context,
												 double *const a_vector,
												 const int64_t size_max);
qmckl_exit_code qmckl_get_jastrow_champ_b_vector(qmckl_context context,
												 double *const b_vector,
												 const int64_t size_max);
qmckl_exit_code qmckl_get_jastrow_champ_c_vector(qmckl_context context,
												 double *const c_vector,
												 const int64_t size_max);
qmckl_exit_code
qmckl_get_jastrow_champ_rescale_factor_ee(const qmckl_context context,
										  double *const rescale_factor_ee);
qmckl_exit_code
qmckl_get_jastrow_champ_rescale_factor_en(const qmckl_context context,
										  double *const rescale_factor_en,
										  const int64_t size_max);
qmckl_exit_code
qmckl_get_jastrow_champ_dim_c_vector(qmckl_context context,
									 int64_t *const dim_c_vector);

//**********
// COMPUTES
//**********

//**********
// PROVIDE
//**********

// Total Jastrow
qmckl_exit_code qmckl_provide_jastrow_champ_value(qmckl_context context);

// Electron/electron component
qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee(qmckl_context context);

// Electron/nucleus component
qmckl_exit_code qmckl_provide_jastrow_champ_factor_en(qmckl_context context);

// Electron/electron/nucleus component
qmckl_exit_code qmckl_provide_jastrow_champ_factor_een(qmckl_context context);

// Distances
qmckl_exit_code qmckl_provide_ee_distance_rescaled(qmckl_context context);
qmckl_exit_code qmckl_provide_en_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_provide_een_rescaled_e(qmckl_context context);
qmckl_exit_code qmckl_provide_een_rescaled_n(qmckl_context context);
qmckl_exit_code
qmckl_provide_jastrow_champ_c_vector_full(qmckl_context context);
qmckl_exit_code qmckl_provide_lkpm_combined_index(qmckl_context context);
qmckl_exit_code qmckl_provide_tmp_c(qmckl_context context);

qmckl_exit_code qmckl_provide_ee_distance(qmckl_context context);
qmckl_exit_code qmckl_provide_en_distance(qmckl_context context);

//**********
// GETTERS (for computes)
//**********

// Total Jastrow
qmckl_exit_code qmckl_get_jastrow_champ_value(qmckl_context context,
											  double *const value,
											  const int64_t size_max);

// Electron/electron component
qmckl_exit_code qmckl_get_jastrow_champ_factor_ee(qmckl_context context,
												  double *const factor_ee,
												  const int64_t size_max);

// Electron/nucleus component
qmckl_exit_code qmckl_get_jastrow_champ_factor_en(qmckl_context context,
												  double *const factor_en,
												  const int64_t size_max);

// Electron/electron/nucleus component
qmckl_exit_code qmckl_get_jastrow_champ_factor_een(qmckl_context context,
												   double *const factor_een,
												   const int64_t size_max);

// Distances
qmckl_exit_code
qmckl_get_jastrow_champ_ee_distance_rescaled(qmckl_context context,
											 double *const distance_rescaled);
qmckl_exit_code
qmckl_get_electron_en_distance_rescaled(qmckl_context context,
										double *distance_rescaled);

qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e(qmckl_context context,
									   double *const distance_rescaled,
									   const int64_t size_max);
qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n(qmckl_context context,
									   double *const distance_rescaled,
									   const int64_t size_max);
qmckl_exit_code qmckl_get_jastrow_champ_tmp_c(qmckl_context context,
											  double *const tmp_c);
