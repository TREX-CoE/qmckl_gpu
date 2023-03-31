#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>
#include "qmckl_basic_type.h"

/* Context */
qmckl_context_device
qmckl_context_check_device(const qmckl_context_device context);

qmckl_exit_code_device
qmckl_context_device_touch_device(const qmckl_context_device context);

void qmckl_lock_device(qmckl_context_device context);
void qmckl_unlock_device(qmckl_context_device context);

/* Error */
qmckl_exit_code_device
qmckl_failwith_device(qmckl_context_device context,
					  const qmckl_exit_code_device exit_code,
					  const char *function, const char *message);

/* Electron */
qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context,
							  const int64_t up_num, const int64_t down_num);

qmckl_exit_code_device qmckl_init_electron_device(qmckl_context_device context);

bool qmckl_nucleus_provided_device(qmckl_context_device context);

/* Point */

qmckl_exit_code_device qmckl_init_point_device(qmckl_context_device context);

/* AO */

qmckl_exit_code_device qmckl_init_ao_basis_device(qmckl_context_device context);
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

/* MO */

qmckl_exit_code_device qmckl_init_mo_basis_device(qmckl_context_device context);
/* Determinant */

qmckl_exit_code_device
qmckl_init_determinant_device(qmckl_context_device context);

/* Jastrow */

qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context_device context);

/* BLAS */

qmckl_vector_device qmckl_vector_alloc_device(qmckl_context_device context,
											  const int64_t size);

qmckl_exit_code_device qmckl_vector_free_device(qmckl_context_device context,
												qmckl_vector_device *vector);

qmckl_matrix_device qmckl_matrix_alloc_device(qmckl_context_device context,
											  const int64_t size1,
											  const int64_t size2);

qmckl_exit_code_device qmckl_matrix_free_device(qmckl_context_device context,
												qmckl_matrix_device *matrix);

qmckl_exit_code_device
qmckl_vector_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_vector_device *vector_out);

qmckl_exit_code_device
qmckl_matrix_of_double_device(const qmckl_context_device context,
							  const double *target, const int64_t size_max,
							  qmckl_matrix_device *matrix);

qmckl_exit_code_device qmckl_transpose_device(qmckl_context_device context,
											  const qmckl_matrix_device A,
											  qmckl_matrix_device At);

qmckl_vector_device
qmckl_vector_of_matrix_device(const qmckl_matrix_device matrix);

qmckl_matrix_device
qmckl_matrix_of_vector_device(const qmckl_vector_device vector,
							  const int64_t size1, const int64_t size2);
