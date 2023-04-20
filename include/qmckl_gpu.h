// This is the header file meant to be included by the users.
// It contains prototypes for all GPU public functions, and definition of
// the _device context variant.

#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/types.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

#ifdef HAVE_CUSPARSE
#include <cuda_runtime_api.h>
#include <cusparse_v2.h>
#endif

//**********
// TYPES
//**********

#define VALID_TAG_DEVICE 0xBEEFFACE
#define INVALID_TAG_DEVICE 0xDEADBEEF

#define QMCKL_DEFAULT_PRECISION_DEVICE 53
#define QMCKL_DEFAULT_RANGE_DEVICE 11

#define QMCKL_MAX_FUN_LEN_DEVICE 256
#define QMCKL_MAX_MSG_LEN_DEVICE 1024

#define QMCKL_TENSOR_ORDER_MAX_DEVICE 16
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

typedef struct qmckl_error_struct_device {

	qmckl_exit_code_device exit_code;
	char function[QMCKL_MAX_FUN_LEN_DEVICE];
	char message[QMCKL_MAX_MSG_LEN_DEVICE];

} qmckl_error_struct_device;

/* Numprec */

typedef struct qmckl_numprec_struct_device {
	uint32_t precision;
	uint32_t range;
} qmckl_numprec_struct_device;

#define QMCKL_NULL_CONTEXT_DEVICE (qmckl_context_device)0

/* BLAS */

typedef struct qmckl_vector_device {
	double *restrict data;
	int64_t size;
} qmckl_vector_device;

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */
typedef struct qmckl_matrix_device {
	double *restrict data;
	int64_t size[2];
} qmckl_matrix_device;

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */
typedef struct qmckl_tensor_device {
	double *restrict data;
	int64_t order;
	int64_t size[QMCKL_TENSOR_ORDER_MAX_DEVICE];
} qmckl_tensor_device;

/* Memory */

typedef struct qmckl_memory_info_struct {
	size_t size;
	void *pointer;
} qmckl_memory_info_struct_device;

typedef struct qmckl_memory_struct_device {
	size_t n_allocated;
	size_t array_size;
	qmckl_memory_info_struct_device *element;
} qmckl_memory_struct_device;

static const qmckl_memory_info_struct_device
	qmckl_memory_info_struct_zero_device = {.size = (size_t)0, .pointer = NULL};

/* Context */

typedef int64_t qmckl_context_device;

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

typedef struct qmckl_jastrow_struct_device {
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
} qmckl_jastrow_struct_device;

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

typedef struct qmckl_context_struct_device {
	/* -- State of the library -- */

	/* Device id, only used w/ OpenMP */
	size_t device_id;

	/* Validity checking */
	uint64_t tag;

	/* Numerical precision */
	qmckl_numprec_struct_device numprec;

	/* Thread lock */
	int lock_count;
	pthread_mutex_t mutex;

	/* Error handling */
	qmckl_error_struct_device error;

	/* Memory allocation */
	qmckl_memory_struct_device memory;
	qmckl_memory_struct_device memory_device;

	/* Current date */
	uint64_t date;

	/* Points */
	qmckl_point_struct_device point;

	/* -- Molecular system -- */
	qmckl_nucleus_struct_device nucleus;
	qmckl_electron_struct_device electron;
	qmckl_ao_basis_struct_device ao_basis;
	qmckl_mo_basis_struct_device mo_basis;
	// TODO When available, add the Jastrow struct
	qmckl_jastrow_struct_device jastrow;
	qmckl_determinant_struct_device det;
	qmckl_local_energy_struct_device local_energy;

	/* To be implemented:
	 */

	/* Pointer to implementation-specific data */

	void *qmckl_extra;

} qmckl_context_struct_device;

//**********
// BASIC FUNCS
//**********

/* Error */
qmckl_exit_code_device
qmckl_failwith_device(qmckl_context_device context,
					  const qmckl_exit_code_device exit_code,
					  const char *function, const char *message);

//**********
// MEMORY
//**********

/* Allocs & frees */
void *qmckl_malloc_host(qmckl_context_device context,
						const qmckl_memory_info_struct_device info);

void *qmckl_malloc_device(qmckl_context_device context, size_t size);

qmckl_exit_code_device qmckl_free_host(qmckl_context_device context,
									   void *const ptr);

qmckl_exit_code_device qmckl_free_device(qmckl_context_device context,
										 void *const ptr);

/* Memcpys */

qmckl_exit_code_device qmckl_memcpy_H2D(qmckl_context_device context,
										void *const dest, void *const src,
										size_t size);
qmckl_exit_code_device qmckl_memcpy_D2H(qmckl_context_device context,
										void *const dest, void *const src,
										size_t size);
qmckl_exit_code_device qmckl_memcpy_D2D(qmckl_context_device context,
										void *const dest, void *const src,
										size_t size);

//**********
// BLAS
//**********

qmckl_vector_device qmckl_vector_alloc_device(qmckl_context_device context,
											  const int64_t size);

qmckl_exit_code_device qmckl_vector_free_device(qmckl_context_device context,
												qmckl_vector_device *vector);

qmckl_exit_code_device
qmckl_vector_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_vector_device *vector_out);

qmckl_matrix_device qmckl_matrix_alloc_device(qmckl_context_device context,
											  const int64_t size1,
											  const int64_t size2);

qmckl_exit_code_device qmckl_matrix_free_device(qmckl_context_device context,
												qmckl_matrix_device *matrix);

qmckl_matrix_device qmckl_matrix_set_device(qmckl_matrix_device matrix,
											double value);

qmckl_exit_code_device
qmckl_matrix_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_matrix_device *matrix_out);

qmckl_exit_code_device qmckl_transpose_device(qmckl_context_device context,
											  const qmckl_matrix_device A,
											  qmckl_matrix_device At);

qmckl_tensor_device qmckl_tensor_alloc_device(qmckl_context_device context,
											  const int64_t order,
											  const int64_t *size);

qmckl_exit_code_device qmckl_tensor_free_device(qmckl_context_device context,
												qmckl_tensor_device *tensor);

qmckl_tensor_device qmckl_tensor_set_device(qmckl_tensor_device tensor,
											double value);

//**********
// CONTEXT
//**********

qmckl_exit_code_device
qmckl_context_touch_device(const qmckl_context_device context);

qmckl_exit_code_device qmckl_init_point_device(qmckl_context_device context);
qmckl_exit_code_device qmckl_init_ao_basis_device(qmckl_context_device context);
qmckl_exit_code_device qmckl_init_mo_basis_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_init_determinant_device(qmckl_context_device context);
qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context_device context);

qmckl_context_device qmckl_context_create_device(int device_id);
qmckl_exit_code_device
qmckl_context_destroy_device(const qmckl_context_device context);

static inline size_t qmckl_get_device_id(qmckl_context_device context) {
	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	return ctx->device_id;
}

qmckl_context_device
qmckl_context_check_device(const qmckl_context_device context);

void qmckl_lock_device(qmckl_context_device context);
void qmckl_unlock_device(qmckl_context_device context);

//**********
// DISTANCE
//**********

qmckl_exit_code_device qmckl_distance_rescaled_device(
	const qmckl_context_device context, const char transa, const char transb,
	const int64_t m, const int64_t n, const double *A, const int64_t lda,
	const double *B, const int64_t ldb, double *const C, const int64_t ldc,
	const double rescale_factor_kappa);

//**********
// POINT
//**********

qmckl_exit_code_device qmckl_set_point_device(qmckl_context_device context,
											  char transp, int64_t num,
											  double *coord, int64_t size_max);

//**********
// NUCLEUS
//**********

bool qmckl_nucleus_provided_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_get_nucleus_num_device(qmckl_context_device context, int64_t *num);
qmckl_exit_code_device
qmckl_get_nucleus_num_device(qmckl_context_device context, int64_t *num);

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num);

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num);

qmckl_exit_code_device
qmckl_set_nucleus_charge_device(qmckl_context_device context, double *charge,
								int64_t size_max);
qmckl_exit_code_device
qmckl_set_nucleus_coord_device(qmckl_context_device context, char transp,
							   double *coord, int64_t size_max);

qmckl_exit_code_device
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_finalize_nucleus_basis_device(qmckl_context_device context);

//**********
// ELECTRON
//**********

bool qmckl_electron_provided_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context, int64_t up_num,
							  int64_t down_num);

qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context, int64_t up_num,
							  int64_t down_num);

qmckl_exit_code_device
qmckl_set_electron_coord_device(qmckl_context_device context, char transp,
								int64_t walk_num, double *coord,
								int64_t size_max);

//**********
// AO
//**********

qmckl_exit_code_device qmckl_compute_ao_basis_shell_gaussian_vgl_device(
	qmckl_context_device context, int prim_num, int shell_num, int point_num,
	int nucl_num, int64_t *nucleus_shell_num, int64_t *nucleus_index,
	double *nucleus_range, int64_t *shell_prim_index, int64_t *shell_prim_num,
	double *coord, double *nucl_coord, double *expo, double *coef_normalized,
	double *shell_vgl);

qmckl_exit_code_device qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_vgl);

qmckl_exit_code_device qmckl_compute_ao_value_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_value);

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_value_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
								 double *const ao_vgl, const int64_t size_max);

qmckl_exit_code_device qmckl_get_ao_basis_ao_value_device(
	qmckl_context_device context, double *const ao_vgl, const int64_t size_max);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_value_inplace_device(qmckl_context_device context,
										   double *const ao_value,
										   const int64_t size_max);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num);
qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(const qmckl_context_device context,
									int64_t *const shell_num);
qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(const qmckl_context_device context,
								   int64_t *const prim_num);
qmckl_exit_code_device
qmckl_get_ao_basis_type_device(const qmckl_context_device context,
							   char *const basis_type);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_finalize_ao_basis_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_set_ao_basis_type_device(qmckl_context_device context, char basis_type);
qmckl_exit_code_device
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num);
qmckl_exit_code_device
qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);
qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_set_ao_basis_ao_num_device(qmckl_context_device context, int64_t ao_num);
qmckl_exit_code_device qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max);
qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max);
qmckl_exit_code_device
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(qmckl_context_device context, int64_t *ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t *ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t *prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_type_device(qmckl_context_device context, char *type);

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_shell_num_device(
	qmckl_context_device context, int64_t *nucleus_shell_num, int64_t nucl_num);
qmckl_exit_code_device qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t nucl_num);

qmckl_exit_code_device qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t shell_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t shell_num);

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t shell_num);

qmckl_exit_code_device
qmckl_get_ao_basis_shell_prim_index_device(qmckl_context_device context,
										   int64_t *nucleus_shell_num,
										   int64_t shell_num);

qmckl_exit_code_device
qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t prim_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t ao_num);

qmckl_exit_code_device
qmckl_get_ao_basis_ao_vgl_inplace_device(qmckl_context_device context,
										 double *ao_vgl, int64_t size_max);

bool qmckl_ao_basis_provided_device(qmckl_context_device context);

//**********
// MO
//**********

bool qmckl_mo_basis_provided_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_num_device(qmckl_context_device context, int64_t *mo_num);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_vgl_device(qmckl_context_device context,
								 double *const mo_vgl, const int64_t size_max);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_value_device(qmckl_context_device context,
								   double *mo_value, int64_t size_max);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_value_inplace_device(qmckl_context_device context,
										   double *mo_value, int64_t size_max);

qmckl_exit_code_device
qmckl_get_mo_basis_mo_vgl_inplace_device(qmckl_context_device context,
										 double *mo_vgl, int64_t size_max);

qmckl_exit_code_device
qmckl_provide_mo_basis_mo_value_device(qmckl_context_device context);

qmckl_exit_code_device
qmckl_provide_mo_basis_mo_vgl_device(qmckl_context_device context);

qmckl_exit_code_device qmckl_compute_mo_basis_mo_value_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t,
	double *restrict ao_value, double *restrict mo_value);

qmckl_exit_code_device qmckl_compute_mo_basis_mo_vgl_device(
	qmckl_context_device context, int64_t ao_num, int64_t mo_num,
	int64_t point_num, double *restrict coefficient_t, double *restrict ao_vgl,
	double *restrict mo_vgl);

qmckl_exit_code_device
qmckl_finalize_mo_basis_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_set_mo_basis_mo_num_device(qmckl_context_device context, int64_t mo_num);
qmckl_exit_code_device
qmckl_set_mo_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient);

//**********
// JASTROW
//**********
// TODO

//****************************
// SHERMAN-MORRISON & WOODBURY
//****************************

// TODO In the future, we could to generate this header on the fly so it
// contains exactly the functions that are enabled (no need for the end user to
// worry about preprocessor)

#ifdef HAVE_CUBLAS
qmckl_exit_code_device qmckl_woodbury_kxk(
	const qmckl_context_device context, cublasHandle_t b_handle,
	cusolverDnHandle_t s_handle, const uint64_t Lds, const uint64_t Dim,
	const uint64_t N_updates, const double *__restrict Updates,
	const uint64_t *__restrict Updates_index, const double breakdown,
	double *__restrict Slater_inv, double *__restrict determinant);
#endif
