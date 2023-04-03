#pragma once

#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>

#include "qmckl_types.h"
#include "qmckl_types.h"

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
