#pragma once

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/* QMCKL DEFINE */

#define VALID_TAG_DEVICE 0xBEEFFACE
#define INVALID_TAG_DEVICE 0xDEADBEEF

#define QMCKL_DEFAULT_PRECISION_DEVICE 53
#define QMCKL_DEFAULT_RANGE_DEVICE 11

#define QMCKL_MAX_FUN_LEN_DEVICE 256
#define QMCKL_MAX_MSG_LEN_DEVICE 1024

#define qmckl_mat_device(m, i, j) m.data[(i) + (j)*m.size[0]]

#define QMCKL_SUCCESS_DEVICE ((qmckl_exit_code_device)0)
#define QMCKL_INVALID_ARG_1_DEVICE ((qmckl_exit_code_device)1)
#define QMCKL_INVALID_ARG_2_DEVICE ((qmckl_exit_code_device)2)
#define QMCKL_INVALID_ARG_3_DEVICE ((qmckl_exit_code_device)3)
#define QMCKL_INVALID_ARG_4_DEVICE ((qmckl_exit_code_device)4)
#define QMCKL_INVALID_ARG_5_DEVICE ((qmckl_exit_code_device)5)
#define QMCKL_INVALID_ARG_6_DEVICE ((qmckl_exit_code_device)6)
#define QMCKL_INVALID_ARG_7_DEVICE ((qmckl_exit_code_device)7)
#define QMCKL_INVALID_ARG_8_DEVICE ((qmckl_exit_code_device)8)
#define QMCKL_INVALID_ARG_9_DEVICE ((qmckl_exit_code_device)9)
#define QMCKL_INVALID_ARG_10_DEVICE ((qmckl_exit_code_device)10)
#define QMCKL_INVALID_ARG_11_DEVICE ((qmckl_exit_code_device)11)
#define QMCKL_INVALID_ARG_12_DEVICE ((qmckl_exit_code_device)12)
#define QMCKL_INVALID_ARG_13_DEVICE ((qmckl_exit_code_device)13)
#define QMCKL_INVALID_ARG_14_DEVICE ((qmckl_exit_code_device)14)
#define QMCKL_INVALID_ARG_15_DEVICE ((qmckl_exit_code_device)15)
#define QMCKL_INVALID_ARG_16_DEVICE ((qmckl_exit_code_device)16)
#define QMCKL_INVALID_ARG_17_DEVICE ((qmckl_exit_code_device)17)
#define QMCKL_INVALID_ARG_18_DEVICE ((qmckl_exit_code_device)18)
#define QMCKL_INVALID_ARG_19_DEVICE ((qmckl_exit_code_device)19)
#define QMCKL_INVALID_ARG_20_DEVICE ((qmckl_exit_code_device)20)
#define QMCKL_FAILURE_DEVICE ((qmckl_exit_code_device)101)
#define QMCKL_ERRNO_DEVICE ((qmckl_exit_code_device)102)
#define QMCKL_INVALID_CONTEXT_DEVICE ((qmckl_exit_code_device)103)
#define QMCKL_ALLOCATION_FAILED_DEVICE ((qmckl_exit_code_device)104)
#define QMCKL_DEALLOCATION_FAILED_DEVICE ((qmckl_exit_code_device)105)
#define QMCKL_NOT_PROVIDED_DEVICE ((qmckl_exit_code_device)106)
#define QMCKL_OUT_OF_BOUNDS_DEVICE ((qmckl_exit_code_device)107)
#define QMCKL_ALREADY_SET_DEVICE ((qmckl_exit_code_device)108)
#define QMCKL_INVALID_EXIT_CODE_DEVICE ((qmckl_exit_code_device)109)


/* Error type */

typedef int32_t qmckl_exit_code_device;

/* #+RESULTS: */
/* :results: */

/* Data structure */

typedef struct qmckl_error_struct_device {

	qmckl_exit_code_device exit_code;
	char function[QMCKL_MAX_FUN_LEN_DEVICE];
	char message[QMCKL_MAX_MSG_LEN_DEVICE];

} qmckl_error_struct_device;

typedef struct qmckl_numprec_struct_device {
	uint32_t precision;
	uint32_t range;
} qmckl_numprec_struct_device;

/* Vector */

/*   | Variable | Type      | Description             | */
/*   |----------+-----------+-------------------------| */
/*   | ~size~   | ~int64_t~ | Dimension of the vector | */
/*   | ~data~   | ~double*~ | Elements                | */

typedef struct qmckl_vector_device {
	double *restrict data;
	int64_t size;
} qmckl_vector_device;

/* Matrix */

/*   | Variable | Type         | Description                 | */
/*   |----------+--------------+-----------------------------| */
/*   | ~size~   | ~int64_t[2]~ | Dimension of each component | */
/*   | ~data~   | ~double*~    | Elements                    | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */

typedef struct qmckl_matrix_device {
	double *restrict data;
	int64_t size[2];
} qmckl_matrix_device;

/* Tensor */

/*   | Variable | Type                              | Description | */
/*   |----------+-----------------------------------+-----------------------------|
 */
/*   | ~order~  | ~int64_t~                         | Order of the tensor | */
/*   | ~size~   | ~int64_t[QMCKL_TENSOR_ORDER_MAX]~ | Dimension of each
 * component | */
/*   | ~data~   | ~double*~                         | Elements | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */

#define QMCKL_TENSOR_ORDER_MAX_DEVICE 16

typedef struct qmckl_local_energy_struct_device {
	double *e_kin;
	double *e_pot;
	double *e_local;
	double *accep_prob;
	double *r_drift;
	double *y_move;
	uint64_t e_kin_date;
	uint64_t e_pot_date;
	uint64_t e_local_date;
	uint64_t accep_prob_date;
	uint64_t r_drift_date;
	uint64_t y_move_date;

	int32_t uninitialized;
	bool provided;
} qmckl_local_energy_struct_device;

typedef struct qmckl_determinant_struct_device {
	char type;
	int64_t det_num_alpha;
	int64_t det_num_beta;
	int64_t up_num;
	int64_t down_num;
	int64_t *mo_index_alpha;
	int64_t *mo_index_beta;

	double *det_value_alpha;
	double *det_value_beta;
	double *det_vgl_alpha;
	double *det_adj_matrix_alpha;
	double *det_inv_matrix_alpha;
	double *det_vgl_beta;
	double *det_adj_matrix_beta;
	double *det_inv_matrix_beta;
	uint64_t det_value_alpha_date;
	uint64_t det_vgl_alpha_date;
	uint64_t det_adj_matrix_alpha_date;
	uint64_t det_inv_matrix_alpha_date;
	uint64_t det_value_beta_date;
	uint64_t det_vgl_beta_date;
	uint64_t det_adj_matrix_beta_date;
	uint64_t det_inv_matrix_beta_date;

	int32_t uninitialized;
	bool provided;
} qmckl_determinant_struct_device;

typedef struct qmckl_jastrow_struct_device_device {
	int32_t uninitialized;
	int64_t aord_num;
	int64_t bord_num;
	int64_t cord_num;
	int64_t type_nucl_num;
	uint64_t asymp_jasa_date;
	uint64_t asymp_jasb_date;
	uint64_t tmp_c_date;
	uint64_t dtmp_c_date;
	uint64_t factor_ee_date;
	uint64_t factor_en_date;
	uint64_t factor_een_date;
	uint64_t factor_ee_deriv_e_date;
	uint64_t factor_en_deriv_e_date;
	uint64_t factor_een_deriv_e_date;
	double rescale_factor_ee;
	double *rescale_factor_en;
	int64_t *type_nucl_vector;
	double *a_vector;
	double *b_vector;
	double *c_vector;
	double *asymp_jasa;
	double *asymp_jasb;
	double *factor_ee;
	double *factor_en;
	double *factor_een;
	double *factor_ee_deriv_e;
	double *factor_en_deriv_e;
	double *factor_een_deriv_e;
	int64_t dim_c_vector;
	uint64_t dim_c_vector_date;
	double *c_vector_full;
	uint64_t c_vector_full_date;
	int64_t *lkpm_combined_index;
	uint64_t lkpm_combined_index_date;
	double *tmp_c;
	double *dtmp_c;
	uint64_t ee_distance_rescaled_date;
	uint64_t ee_distance_rescaled_deriv_e_date;
	uint64_t en_distance_rescaled_date;
	uint64_t en_distance_rescaled_deriv_e_date;
	double *ee_distance_rescaled;
	double *ee_distance_rescaled_deriv_e;
	double *en_distance_rescaled;
	double *en_distance_rescaled_deriv_e;
	double *een_rescaled_e;
	double *een_rescaled_n;
	uint64_t een_rescaled_e_date;
	uint64_t een_rescaled_n_date;
	double *een_rescaled_e_deriv_e;
	double *een_rescaled_n_deriv_e;
	uint64_t een_rescaled_e_deriv_e_date;
	uint64_t een_rescaled_n_deriv_e_date;
	bool provided;
	char *type;

	bool gpu_offload;
} qmckl_jastrow_struct;

typedef struct qmckl_mo_basis_struct_device {
	int64_t mo_num;
	double *restrict coefficient;
	double *restrict coefficient_t;

	double *restrict mo_vgl;
	double *restrict mo_value;
	uint64_t mo_vgl_date;
	uint64_t mo_value_date;

	int32_t uninitialized;
	bool provided;
} qmckl_mo_basis_struct_device;

typedef struct qmckl_tensor_device {
	double *restrict data;
	int64_t order;
	int64_t size[QMCKL_TENSOR_ORDER_MAX_DEVICE];
} qmckl_tensor_device;

typedef struct qmckl_nucleus_struct_device {
	int64_t num;
	int64_t repulsion_date;
	int64_t nn_distance_date;
	int64_t coord_date;
	qmckl_vector_device charge;
	qmckl_matrix_device coord;
	qmckl_matrix_device nn_distance;
	double repulsion;
	int32_t uninitialized;
	bool provided;
} qmckl_nucleus_struct_device;

typedef struct qmckl_ao_basis_struct_device {
	int64_t shell_num;
	int64_t prim_num;
	int64_t ao_num;
	int64_t *restrict nucleus_index;
	int64_t *restrict nucleus_shell_num;
	int32_t *restrict shell_ang_mom;
	int64_t *restrict shell_prim_num;
	int64_t *restrict shell_prim_index;
	double *restrict shell_factor;
	double *restrict exponent;
	double *restrict coefficient;
	double *restrict prim_factor;
	double *restrict ao_factor;

	int64_t *restrict nucleus_prim_index;
	double *restrict coefficient_normalized;
	int32_t *restrict nucleus_max_ang_mom;
	double *restrict nucleus_range;
	double *restrict primitive_vgl;
	uint64_t primitive_vgl_date;
	double *restrict shell_vgl;
	uint64_t shell_vgl_date;
	double *restrict ao_vgl;
	uint64_t ao_vgl_date;
	double *restrict ao_value;
	uint64_t ao_value_date;

	int32_t uninitialized;
	bool provided;
	bool ao_cartesian;
	char type;
	/* HPC specific data structures */
	int32_t *restrict prim_num_per_nucleus;
	qmckl_tensor_device coef_per_nucleus;
	qmckl_matrix_device expo_per_nucleus;
} qmckl_ao_basis_struct_device;

typedef struct qmckl_point_struct_device {
	int64_t num;
	uint64_t date;
	qmckl_matrix_device coord;
} qmckl_point_struct_device;

typedef int64_t qmckl_context_device;
#define QMCKL_NULL_CONTEXT_DEVICE (qmckl_context_device)0

typedef struct qmckl_walker_struct_device {
	int64_t num;
	qmckl_point_struct_device point;
} qmckl_walker_device;

typedef struct qmckl_electron_struct_device {
	int64_t num;
	int64_t up_num;
	int64_t down_num;
	qmckl_walker_device walker;
	qmckl_walker_device walker_old;
	uint64_t ee_distance_date;
	uint64_t en_distance_date;
	uint64_t ee_potential_date;
	uint64_t en_potential_date;
	double *ee_distance;
	double *en_distance;
	double *ee_potential;
	double *en_potential;
	int32_t uninitialized;
	bool provided;
} qmckl_electron_struct_device;
