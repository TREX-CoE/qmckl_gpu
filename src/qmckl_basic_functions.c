#include "include/qmckl_basic_functions.h"

#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <errno.h>
#include <math.h>

/* Error */

qmckl_exit_code_device qmckl_set_error_device(qmckl_context_device context,
								const qmckl_exit_code_device exit_code,
								const char *function_name,
								const char *message) {
	/* Passing a function name and a message is mandatory. */
	assert(function_name != NULL);
	assert(message != NULL);

	/* Exit codes are assumed valid. */
	assert(exit_code >= 0);
	assert(exit_code != QMCKL_SUCCESS);
	assert(exit_code < QMCKL_INVALID_EXIT_CODE);

	/* The context is assumed to exist. */
	assert(qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

	qmckl_lock_device(context);
	{
		qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
		assert(ctx != NULL); /* Impossible because the context is valid. */

		ctx->error.exit_code = exit_code;
		strncpy(ctx->error.function, function_name, QMCKL_MAX_FUN_LEN - 1);
		strncpy(ctx->error.message, message, QMCKL_MAX_MSG_LEN - 1);
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS;
}

const char *qmckl_string_of_error_device(const qmckl_exit_code_device error) {
	switch (error) {}
	return "Unknown error";
}

qmckl_exit_code_device qmckl_failwith_device(qmckl_context_device context,
							   const qmckl_exit_code_device exit_code,
							   const char *function, const char *message) {
	assert(exit_code > 0);
	assert(exit_code < QMCKL_INVALID_EXIT_CODE);
	assert(function != NULL);
	assert(strlen(function) < QMCKL_MAX_FUN_LEN);
	if (message != NULL) {
		assert(strlen(message) < QMCKL_MAX_MSG_LEN);
	}

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
		return QMCKL_INVALID_CONTEXT;

	if (message == NULL) {
		qmckl_exit_code_device rc = qmckl_set_error(context, exit_code, function,
											 qmckl_string_of_error(exit_code));
		assert(rc == QMCKL_SUCCESS);
	} else {
		qmckl_exit_code_device rc =
			qmckl_set_error(context, exit_code, function, message);
		assert(rc == QMCKL_SUCCESS);
	}

	return exit_code;
}

/* Electron */

qmckl_exit_code_device qmckl_init_electron_device(qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->electron.uninitialized = (1 << 1) - 1;

	return QMCKL_SUCCESS;
}

bool qmckl_electron_provided_device(const qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	return ctx->electron.provided;
}

qmckl_exit_code_device qmckl_set_electron_num_device(qmckl_context_device context,
									   const int64_t up_num,
									   const int64_t down_num) {
	int32_t mask = 1 << 0;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET,
							  "qmckl_set_electron_*", NULL);
	}

	if (up_num <= 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_electron_num", "up_num <= 0");
	}

	if (down_num < 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_electron_num", "down_num < 0");
	}

	ctx->electron.up_num = up_num;
	ctx->electron.down_num = down_num;
	ctx->electron.num = up_num + down_num;

	ctx->electron.uninitialized &= ~mask;
	ctx->electron.provided = (ctx->electron.uninitialized == 0);

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_set_electron_coord_device(qmckl_context_device context,
										 const char transp,
										 const int64_t walk_num,
										 const double *coord,
										 const int64_t size_max) {

	int32_t mask = 0;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET,
							  "qmckl_set_electron_*", NULL);
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_electron_coord",
							  "transp should be 'N' or 'T'");
	}

	if (walk_num <= 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_electron_coord", "walk_num <= 0");
	}

	if (coord == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_set_electron_coord",
							  "coord is a null pointer");
	}

	const int64_t elec_num = ctx->electron.num;

	if (elec_num == 0L) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_electron_coord",
							  "elec_num is not set");
	}

	/* Swap pointers */
	qmckl_walker tmp = ctx->electron.walker_old;
	ctx->electron.walker_old = ctx->electron.walker;
	ctx->electron.walker = tmp;

	memcpy(&(ctx->point), &(ctx->electron.walker.point),
		   sizeof(qmckl_point_struct));
	qmckl_exit_code_device rc;
	rc = qmckl_set_point(context, transp, walk_num * elec_num, coord, size_max);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->electron.walker.num = walk_num;
	memcpy(&(ctx->electron.walker.point), &(ctx->point),
		   sizeof(qmckl_point_struct));

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_electron_num_device(const qmckl_context_device context,
									   int64_t *const num) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (num == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_electron_num",
							  "num is a null pointer");
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->electron.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED;
	}

	assert(ctx->electron.num > (int64_t)0);
	*num = ctx->electron.num;
	return QMCKL_SUCCESS;
}

/* nucleus */

qmckl_exit_code_device qmckl_get_nucleus_num_device(const qmckl_context_device context,
									  int64_t *const num) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (num == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_nucleus_num", "num is a null pointer");
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->nucleus.uninitialized & mask) != 0) {
		*num = (int64_t)0;
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_nucleus_num",
							  "nucleus data is not provided");
	}

	assert(ctx->nucleus.num >= (int64_t)0);
	*num = ctx->nucleus.num;

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_init_nucleus_device(qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->nucleus.uninitialized = (1 << 3) - 1;

	/* Default values */

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_set_nucleus_num_device(qmckl_context_device context,
									  const int64_t num) {
	int32_t mask = 1 << 0;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET, "qmckl_set_nucleus_*",
							  NULL);
	}

	if (num <= 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_nucleus_num", "num <= 0");
	}

	ctx->nucleus.num = num;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

/* Sets the nuclear charges of all the atoms. */

qmckl_exit_code_device qmckl_set_nucleus_charge_device_device(qmckl_context_device context,
										 const double *charge,
										 const int64_t size_max) {
	int32_t mask = 1 << 1;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET, "qmckl_set_nucleus_*",
							  NULL);
	}

	if (charge == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_nucleus_charge",
							  "charge is a null pointer");
	}

	int64_t num;
	qmckl_exit_code_device rc;

	rc = qmckl_get_nucleus_num(context, &num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	if (num > size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_nucleus_charge", "Array too small");
	}

	ctx->nucleus.charge = qmckl_vector_alloc(context, num);
	rc = qmckl_vector_of_double(context, charge, num, &(ctx->nucleus.charge));
	if (rc != QMCKL_SUCCESS) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_nucleus_charge",
							  "Error in vector->double* conversion");
	}

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_set_nucleus_coord_device(qmckl_context_device context,
										const char transp, const double *coord,
										const int64_t size_max) {
	int32_t mask = 1 << 2;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET, "qmckl_set_nucleus_*",
							  NULL);
	}

	qmckl_exit_code_device rc;

	const int64_t nucl_num = (int64_t)ctx->nucleus.num;

	if (ctx->nucleus.coord.data != NULL) {
		rc = qmckl_matrix_free(context, &(ctx->nucleus.coord));
		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	ctx->nucleus.coord = qmckl_matrix_alloc(context, nucl_num, 3);

	if (ctx->nucleus.coord.data == NULL) {
		return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_nucleus_coord", NULL);
	}

	if (size_max < 3 * nucl_num) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_set_nucleus_coord", "Array too small");
	}

	if (transp == 'N') {
		qmckl_matrix At;
		At = qmckl_matrix_alloc(context, 3, nucl_num);
		rc = qmckl_matrix_of_double(context, coord, 3 * nucl_num, &At);
		if (rc != QMCKL_SUCCESS)
			return rc;
		rc = qmckl_transpose(context, At, ctx->nucleus.coord);
	} else {
		rc = qmckl_matrix_of_double(context, coord, nucl_num * 3,
									&(ctx->nucleus.coord));
	}
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

bool qmckl_nucleus_provided_device(qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	return ctx->nucleus.provided;
}

/* Memory */

qmckl_exit_code_device qmckl_get_malloc_info_device(qmckl_context context, const void *ptr,
									  qmckl_memory_info_struct *info) {

	assert(qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (ptr == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_get_malloc_info", "Null pointer");
	}

	if (info == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_malloc_info", "Null pointer");
	}

	qmckl_lock(context);
	{
		/* Find the pointer entry */
		size_t pos = (size_t)0;
		while (pos < ctx->memory.array_size &&
			   ctx->memory.element[pos].pointer != ptr) {
			pos += (size_t)1;
		}

		if (pos >= ctx->memory.array_size) {
			/* Not found */
			qmckl_unlock(context);
			return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
								  "qmckl_get_malloc_info",
								  "Pointer not found in context");
		}

		/* Copy info */
		memcpy(info, &(ctx->memory.element[pos]),
			   sizeof(qmckl_memory_info_struct));
	}
	qmckl_unlock(context);

	return QMCKL_SUCCESS;
}

/* Point */

qmckl_exit_code_device qmckl_set_point_device(qmckl_context context, const char transp,
								const int64_t num, const double *coord,
								const int64_t size_max) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	if (num <= 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_set_point",
							  "Number of points should be >0.");
	}

	if (size_max < 3 * num) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4, "qmckl_set_point",
							  "Array too small");
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_set_point",
							  "transp should be 'N' or 'T'");
	}

	if (coord == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_set_point",
							  "coord is a NULL pointer");
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	qmckl_exit_code_device rc;
	if (num != ctx->point.num) {

		if (ctx->point.coord.data != NULL) {
			rc = qmckl_matrix_free(context, &(ctx->point.coord));
			assert(rc == QMCKL_SUCCESS);
		}

		ctx->point.coord = qmckl_matrix_alloc(context, num, 3);
		if (ctx->point.coord.data == NULL) {
			return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
								  "qmckl_set_point", NULL);
		}
	};

	ctx->point.num = num;
	if (transp == 'T') {
		double *a = ctx->point.coord.data;
#ifdef HAVE_OPENMP
#pragma omp for
#endif
		for (int64_t i = 0; i < 3 * num; ++i) {
			a[i] = coord[i];
		}
	} else {
#ifdef HAVE_OPENMP
#pragma omp for
#endif
		for (int64_t i = 0; i < num; ++i) {
			qmckl_mat(ctx->point.coord, i, 0) = coord[3 * i];
			qmckl_mat(ctx->point.coord, i, 1) = coord[3 * i + 1];
			qmckl_mat(ctx->point.coord, i, 2) = coord[3 * i + 2];
		}
	}

	/* Increment the date of the context */
	rc = qmckl_context_touch(context);
	assert(rc == QMCKL_SUCCESS);

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_init_point_device(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	memset(&(ctx->point), 0, sizeof(qmckl_point_struct));

	return QMCKL_SUCCESS;
}

/* AO */

qmckl_exit_code_device qmckl_init_ao_basis_device(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_init_ao_basis", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->ao_basis.uninitialized = (1 << 14) - 1;

	/* Default values */
	ctx->ao_basis.ao_cartesian = true;

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_ao_basis_ao_num_device(const qmckl_context context,
										  int64_t *const ao_num) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_num", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 12;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_ao_num", NULL);
	}

	assert(ctx->ao_basis.ao_num > (int64_t)0);

	*ao_num = ctx->ao_basis.ao_num;
	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_num_device(const qmckl_context context,
											 int64_t *const shell_num) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_factor", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_shell_num", NULL);
	}

	assert(ctx->ao_basis.shell_num > (int64_t)0);
	*shell_num = ctx->ao_basis.shell_num;
	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_ao_basis_prim_num_device(const qmckl_context context,
											int64_t *const prim_num) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_prim_num", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 2;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_prim_num", NULL);
	}

	assert(ctx->ao_basis.prim_num > (int64_t)0);

	*prim_num = ctx->ao_basis.prim_num;
	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_ao_basis_type_device(const qmckl_context context,
										char *const basis_type) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_type", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_ao_basis_type", NULL);
	}

	assert(ctx->ao_basis.type != (char)0);

	basis_type[0] = ctx->ao_basis.type;
	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_get_mo_basis_mo_num_device(const qmckl_context context,
										  int64_t *mo_num) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_mo_basis_mo_num", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->mo_basis.uninitialized & mask) != 0) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_get_mo_basis_mo_num", NULL);
	}

	assert(ctx->mo_basis.mo_num > (int64_t)0);
	*mo_num = ctx->mo_basis.mo_num;
	return QMCKL_SUCCESS;
}

/* MO */

qmckl_exit_code_device qmckl_init_mo_basis_device(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->mo_basis.uninitialized = (1 << 2) - 1;

	return QMCKL_SUCCESS;
}

/* Determinant */

qmckl_exit_code_device qmckl_init_determinant_device(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->det.uninitialized = (1 << 5) - 1;

	return QMCKL_SUCCESS;
}

/* Jastrow */

qmckl_exit_code_device qmckl_init_jastrow_device(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	ctx->jastrow.uninitialized = (1 << 10) - 1;

	/* Default values */
	ctx->jastrow.aord_num = -1;
	ctx->jastrow.bord_num = -1;
	ctx->jastrow.cord_num = -1;
	ctx->jastrow.type_nucl_num = -1;
	ctx->jastrow.dim_c_vector = -1;

	return QMCKL_SUCCESS;
}

/* Vector */

qmckl_vector qmckl_vector_alloc_device(qmckl_context context, const int64_t size) {
	/* Should always be true by contruction */
	assert(size > (int64_t)0);

	qmckl_vector result;
	result.size = size;

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = size * sizeof(double);
	result.data = (double *)qmckl_malloc(context, mem_info);

	if (result.data == NULL) {
		result.size = (int64_t)0;
	}

	return result;
}

qmckl_exit_code_device qmckl_vector_free_device(qmckl_context context, qmckl_vector *vector) {
	if (vector == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_vector_free",
							  "Null pointer");
	}

	/* Always true */
	assert(vector->data != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_free(context, vector->data);
	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	vector->size = (int64_t)0;
	vector->data = NULL;
	return QMCKL_SUCCESS;
}

qmckl_matrix qmckl_matrix_alloc_device(qmckl_context context, const int64_t size1,
								const int64_t size2) {
	/* Should always be true by contruction */
	assert(size1 * size2 > (int64_t)0);

	qmckl_matrix result;

	result.size[0] = size1;
	result.size[1] = size2;

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = size1 * size2 * sizeof(double);
	result.data = (double *)qmckl_malloc(context, mem_info);

	if (result.data == NULL) {
		result.size[0] = (int64_t)0;
		result.size[1] = (int64_t)0;
	}

	return result;
}

qmckl_exit_code_device qmckl_matrix_free_device(qmckl_context context, qmckl_matrix *matrix) {
	if (matrix == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_matrix_free",
							  "Null pointer");
	}

	/* Always true */
	assert(matrix->data != NULL);

	qmckl_exit_code_device rc;

	rc = qmckl_free(context, matrix->data);
	if (rc != QMCKL_SUCCESS) {
		return rc;
	}
	matrix->data = NULL;
	matrix->size[0] = (int64_t)0;
	matrix->size[1] = (int64_t)0;

	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_vector_of_double_device(const qmckl_context context,
									   const double *target,
									   const int64_t size_max,
									   qmckl_vector *vector_out) {
	qmckl_vector vector = *vector_out;
	/* Always true by construction */
	assert(qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

	if (vector.size == 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_double_of_vector", "Vector not allocated");
	}

	if (vector.size != size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_double_of_vector", "Wrong vector size");
	}

	for (int64_t i = 0; i < vector.size; ++i) {
		vector.data[i] = target[i];
	}

	*vector_out = vector;
	return QMCKL_SUCCESS;
}

qmckl_exit_code_device qmckl_matrix_of_double_device(const qmckl_context context,
									   const double *target,
									   const int64_t size_max,
									   qmckl_matrix *matrix) {
	qmckl_vector vector = qmckl_vector_of_matrix(*matrix);
	qmckl_exit_code_device rc =
		qmckl_vector_of_double(context, target, size_max, &vector);
	*matrix = qmckl_matrix_of_vector(vector, matrix->size[0], matrix->size[1]);
	return rc;
}

qmckl_exit_code_device qmckl_transpose_device(qmckl_context context, const qmckl_matrix A,
								qmckl_matrix At) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (A.size[0] < 1) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_transpose",
							  "Invalid size for A");
	}

	if (At.data == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_transpose",
							  "Output matrix not allocated");
	}

	if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_transpose",
							  "Invalid size for At");
	}

	for (int64_t j = 0; j < At.size[1]; ++j)
		for (int64_t i = 0; i < At.size[0]; ++i)
			qmckl_mat(At, i, j) = qmckl_mat(A, j, i);

	return QMCKL_SUCCESS;
}

qmckl_vector qmckl_vector_of_matrix_device(const qmckl_matrix matrix) {
	qmckl_vector result;

	result.size = matrix.size[0] * matrix.size[1];
	result.data = matrix.data;

	return result;
}

qmckl_matrix qmckl_matrix_of_vector_device(const qmckl_vector vector,
									const int64_t size1, const int64_t size2) {
	/* Always true */
	assert(size1 * size2 == vector.size);

	qmckl_matrix result;

	result.size[0] = size1;
	result.size[1] = size2;
	result.data = vector.data;

	return result;
}

/* Trexio */

trexio_t *qmckl_trexio_open_X_device(const char *file_name, qmckl_exit_code_device *rc) {
	*rc = QMCKL_SUCCESS;
	trexio_t *file = NULL;

	file = trexio_open(file_name, 'r', TREXIO_TEXT, rc);
	if (file != NULL)
		return file;

	file = trexio_open(file_name, 'r', TREXIO_HDF5, rc);
	if (file != NULL)
		return file;

	*rc = QMCKL_FAILURE;
	/* TODO
	  return qmckl_failwith( context,
							 QMCKL_FAILURE,
							 "trexio_read_electron_up_num",
							 trexio_string_of_error(rcio));
							 */
	return NULL;
}
