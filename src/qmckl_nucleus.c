#include "../include/qmckl_nucleus.h"

qmckl_exit_code_device
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context) {

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;

	//int device_id = qmckl_get_device_id(context);

	ctx->ao_basis.prim_num_per_nucleus = (int32_t *)qmckl_malloc_device(
		context, ctx->nucleus.num * sizeof(int32_t));

	/* Find max number of primitives per nucleus */

	// Extract arrays from context
	int64_t *nucleus_shell_num = ctx->ao_basis.nucleus_shell_num;
	int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
	int64_t *shell_prim_num = ctx->ao_basis.shell_prim_num;
	int32_t *prim_num_per_nucleus = ctx->ao_basis.prim_num_per_nucleus;

	int64_t shell_max = 0;
	int64_t prim_max = 0;
	int64_t nucl_num = ctx->nucleus.num;

	int64_t *shell_max_ptr = &shell_max;
	int64_t *prim_max_ptr = &prim_max;

	// TODO Manually specify OpenACC clauses

#define DEV_PTRS_K01 nucleus_shell_num, nucleus_index, \
                     shell_prim_num, prim_num_per_nucleus

#pragma acc kernels \
            copy(shell_max_ptr[:1], prim_max_ptr[:1]) \
			deviceptr( DEV_PTRS_K01 ) 
#pragma omp target \
            map(tofrom : shell_max_ptr[:1], prim_max_ptr[:1])    \
	        is_device_ptr( DEV_PTRS_K01 )
{

	for (int inucl = 0; inucl < nucl_num; ++inucl) {
		shell_max_ptr[0] = nucleus_shell_num[inucl] > shell_max_ptr[0]
							   ? nucleus_shell_num[inucl]
							   : shell_max_ptr[0];

		int64_t prim_num = 0;
		const int64_t ishell_start = nucleus_index[inucl];
		const int64_t ishell_end =
			nucleus_index[inucl] + nucleus_shell_num[inucl];
		for (int64_t ishell = ishell_start; ishell < ishell_end;
			 ++ishell) {
			prim_num += shell_prim_num[ishell];
		}
		prim_max_ptr[0] =
			prim_num > prim_max_ptr[0] ? prim_num : prim_max_ptr[0];
		prim_num_per_nucleus[inucl] = prim_num;
	}

}

	int64_t size[3] = {prim_max, shell_max, nucl_num};
	ctx->ao_basis.coef_per_nucleus =
		qmckl_tensor_alloc_device(context, 3, size);
	ctx->ao_basis.coef_per_nucleus =
		qmckl_tensor_set_device(ctx->ao_basis.coef_per_nucleus, 0.);

	ctx->ao_basis.expo_per_nucleus =
		qmckl_matrix_alloc_device(context, prim_max, nucl_num);
	ctx->ao_basis.expo_per_nucleus =
		qmckl_matrix_set_device(ctx->ao_basis.expo_per_nucleus, 0.);

	// To avoid offloading structures, expo is split in two arrays :
	// struct combined expo[prim_max];
	// ... gets replaced by :
	double *expo_expo = qmckl_malloc_device(context, prim_max * sizeof(double));
	int64_t *expo_index =
		qmckl_malloc_device(context, prim_max * sizeof(double));

	double *coef =
		qmckl_malloc_device(context, shell_max * prim_max * sizeof(double));
	double *newcoef = qmckl_malloc_device(context, prim_max * sizeof(double));

	int64_t *newidx = qmckl_malloc_device(context, prim_max * sizeof(int64_t));

	int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;
	double *exponent = ctx->ao_basis.exponent;
	double *coefficient_normalized = ctx->ao_basis.coefficient_normalized;

	double *expo_per_nucleus_data = ctx->ao_basis.expo_per_nucleus.data;
	int expo_per_nucleus_s0 = ctx->ao_basis.expo_per_nucleus.size[0];

	double *coef_per_nucleus_data = ctx->ao_basis.coef_per_nucleus.data;
	int coef_per_nucleus_s0 = ctx->ao_basis.coef_per_nucleus.size[0];
	int coef_per_nucleus_s1 = ctx->ao_basis.coef_per_nucleus.size[1];

	// TODO Manually specify OpenACC clauses

#define DEV_PTRS_K02 expo_expo, expo_index, coef, newcoef,            \
                     nucleus_index, shell_prim_index, nucleus_shell_num, \
                     exponent, coefficient_normalized, shell_prim_num,   \
                     expo_per_nucleus_data, coef_per_nucleus_data,       \
                     prim_num_per_nucleus, newidx

#pragma acc kernels deviceptr( DEV_PTRS_K02 )
#pragma omp target is_device_ptr( DEV_PTRS_K02 )
{
	for (int64_t inucl = 0; inucl < nucl_num; ++inucl) {
		for (int i = 0; i < prim_max; i++) {
			expo_expo[i] = 0.;
			expo_index[i] = 0;
		}
		for (int i = 0; i < shell_max * prim_max; i++) {
			coef[i] = 0.;
		}

		int64_t idx = 0;
		const int64_t ishell_start = nucleus_index[inucl];
		const int64_t ishell_end =
			nucleus_index[inucl] + nucleus_shell_num[inucl];

		for (int64_t ishell = ishell_start; ishell < ishell_end;
			 ++ishell) {

			const int64_t iprim_start = shell_prim_index[ishell];
			const int64_t iprim_end =
				shell_prim_index[ishell] + shell_prim_num[ishell];
			for (int64_t iprim = iprim_start; iprim < iprim_end;
				 ++iprim) {
				expo_expo[idx] = exponent[iprim];
				expo_index[idx] = idx;
				idx += 1;
			}
		}

		/* Sort exponents */
		// In the CPU version :
		// qsort( expo, (size_t) idx, sizeof(struct combined),
		// compare_basis );
		// ... is replaced by a hand written bubble sort on
		// expo_expo :
		double tmp;
		for (int i = 0; i < idx - 1; i++) {
			for (int j = 0; j < idx - i - 1; j++) {
				if (expo_expo[j + 1] < expo_expo[j]) {
					tmp = expo_expo[j + 1];
					expo_expo[j + 1] = expo_expo[j];
					expo_expo[j] = tmp;

					tmp = expo_index[j + 1];
					expo_index[j + 1] = expo_index[j];
					expo_index[j] = tmp;
				}
			}
		}

		idx = 0;
		int64_t idx2 = 0;
		for (int64_t ishell = ishell_start; ishell < ishell_end;
			 ++ishell) {

			for (int i = 0; i < prim_max; i++) {
				newcoef[i] = 0;
			}
			const int64_t iprim_start = shell_prim_index[ishell];
			const int64_t iprim_end =
				shell_prim_index[ishell] + shell_prim_num[ishell];

			for (int64_t iprim = iprim_start; iprim < iprim_end;
				 ++iprim) {
				newcoef[idx] = coefficient_normalized[iprim];
				idx += 1;
			}
			for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
				idx2 = expo_index[i];
				coef[(ishell - ishell_start) * prim_max + i] =
					newcoef[idx2];
			}
		}

		/* Apply ordering to coefficients */

		/* Remove duplicates */
		int64_t idxmax = 0;
		idx = 0;
		newidx[0] = 0;

		for (int32_t i = 1; i < prim_num_per_nucleus[inucl]; ++i) {
			if (expo_expo[i] != expo_expo[i - 1]) {
				idx += 1;
			}
			newidx[i] = idx;
		}
		idxmax = idx;

		for (int32_t j = 0; j < ishell_end - ishell_start; ++j) {
			for (int i = 0; i < prim_max; i++) {
				newcoef[i] = 0.;
			}

			for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
				newcoef[newidx[i]] += coef[j * prim_max + i];
			}
			for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
				coef[j * prim_max + i] = newcoef[i];
			}
		}

		for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
			expo_expo[newidx[i]] = expo_expo[i];
		}
		prim_num_per_nucleus[inucl] = (int32_t)idxmax + 1;

		for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
			expo_per_nucleus_data[i + inucl * expo_per_nucleus_s0] =
				expo_expo[i];
		}

		for (int32_t j = 0; j < ishell_end - ishell_start; ++j) {
			for (int32_t i = 0; i < prim_num_per_nucleus[inucl]; ++i) {
				coef_per_nucleus_data[(i) +
									  coef_per_nucleus_s0 *
										  ((j) + coef_per_nucleus_s1 *
													 (inucl))] =
					coef[j * prim_max + i];
			}
		}
	}
	
}
// End of target region

	qmckl_free_device(context, expo_expo);
	qmckl_free_device(context, expo_index);
	qmckl_free_device(context, coef);
	qmckl_free_device(context, newcoef);
	qmckl_free_device(context, newidx);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_finalize_nucleus_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_finalize_nucleus_basis_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t nucl_num = 0;

	qmckl_exit_code_device rc =
		qmckl_get_nucleus_num_device(context, &nucl_num);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	/* nucleus_prim_index */
	{
		ctx->ao_basis.nucleus_prim_index = (int64_t *)qmckl_malloc_device(
			context, (ctx->nucleus.num + (int64_t)1) * sizeof(int64_t));

		if (ctx->ao_basis.nucleus_prim_index == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
										 "ao_basis.nucleus_prim_index", NULL);
		}

		// Extract arrays from context
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_prim_index = ctx->ao_basis.nucleus_prim_index;
		int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;

		int prim_num = ctx->ao_basis.prim_num;

		// TODO Manually specify OpenACC clauses
		
#define DEV_PTRS_K03 nucleus_index, nucleus_prim_index, shell_prim_index

#pragma acc kernels deviceptr( DEV_PTRS_K03 )
#pragma omp target is_device_ptr( DEV_PTRS_K03 )
{
			#pragma omp parallel for
			for (int64_t i = 0; i < nucl_num; ++i) {
				int64_t shell_idx = nucleus_index[i];
				nucleus_prim_index[i] = shell_prim_index[shell_idx];
			}

			nucleus_prim_index[nucl_num] = prim_num;
		
}
	}

	/* Normalize coefficients */
	{

		ctx->ao_basis.coefficient_normalized = (double *)qmckl_malloc_device(
			context, ctx->ao_basis.prim_num * sizeof(double));

		if (ctx->ao_basis.coefficient_normalized == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
				"ao_basis.coefficient_normalized", NULL);
		}

		// Extract arrays from context
		int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;
		int64_t *shell_prim_num = ctx->ao_basis.shell_prim_num;
		double *coefficient_normalized = ctx->ao_basis.coefficient_normalized;
		double *coefficient = ctx->ao_basis.coefficient;
		double *prim_factor = ctx->ao_basis.prim_factor;
		double *shell_factor = ctx->ao_basis.shell_factor;

		int shell_num = ctx->ao_basis.shell_num;

#define DEV_PTRS_K04 shell_prim_index, shell_prim_num,                   \
                     coefficient_normalized, coefficient, prim_factor,   \
                     shell_factor

#pragma acc kernels deviceptr( DEV_PTRS_K04 )
#pragma omp target is_device_ptr( DEV_PTRS_K04 )
{
		for (int64_t ishell = 0; ishell < shell_num; ++ishell) {
			for (int64_t iprim = shell_prim_index[ishell];
				 iprim <
				 shell_prim_index[ishell] + shell_prim_num[ishell];
				 ++iprim) {
				coefficient_normalized[iprim] = coefficient[iprim] *
												prim_factor[iprim] *
												shell_factor[ishell];
			}
		}
}
	}

	/* Find max angular momentum on each nucleus */
	{
		ctx->ao_basis.nucleus_max_ang_mom = (int32_t *)qmckl_malloc_device(
			context, ctx->nucleus.num * sizeof(int32_t));

		if (ctx->ao_basis.nucleus_max_ang_mom == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
										 "ao_basis.nucleus_max_ang_mom", NULL);
		}

		// Extract arrays from context
		int32_t *nucleus_max_ang_mom = ctx->ao_basis.nucleus_max_ang_mom;
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_shell_num = ctx->ao_basis.nucleus_shell_num;
		int32_t *shell_ang_mom = ctx->ao_basis.shell_ang_mom;

#define DEV_PTRS_K05 nucleus_max_ang_mom, nucleus_index, \
                     nucleus_shell_num, shell_ang_mom

#pragma acc kernels deviceptr( DEV_PTRS_K05 ) 
#pragma omp target is_device_ptr( DEV_PTRS_K05 )
{
		#pragma omp parallel for
		for (int64_t inucl = 0; inucl < nucl_num; ++inucl) {
			nucleus_max_ang_mom[inucl] = 0;
			for (int64_t ishell = nucleus_index[inucl];
				 ishell <
				 nucleus_index[inucl] + nucleus_shell_num[inucl];
				 ++ishell) {
				nucleus_max_ang_mom[inucl] =
					nucleus_max_ang_mom[inucl] > shell_ang_mom[ishell]
						? nucleus_max_ang_mom[inucl]
						: shell_ang_mom[ishell];
			}
		}
}
	}

	/* Find distance beyond which all AOs are zero.
	   The distance is obtained by sqrt(log(cutoff)*range) */
	{
		if (ctx->ao_basis.type == 'G') {
			ctx->ao_basis.nucleus_range = (double *)qmckl_malloc_device(
				context, ctx->nucleus.num * sizeof(double));

			if (ctx->ao_basis.nucleus_range == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
											 "ao_basis.nucleus_range", NULL);
			}

			// Extract arrays from context
			double *nucleus_range = ctx->ao_basis.nucleus_range;
			int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
			int64_t *nucleus_shell_num = ctx->ao_basis.nucleus_shell_num;
			int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;
			int64_t *shell_prim_num = ctx->ao_basis.shell_prim_num;
			double *exponent = ctx->ao_basis.exponent;

			int nucleus_num = ctx->nucleus.num;

#define DEV_PTRS_K06 nucleus_range, nucleus_index,        \
                     nucleus_shell_num, shell_prim_index, \
                     shell_prim_num, exponent

#pragma acc kernels deviceptr( DEV_PTRS_K06 ) 
#pragma omp target is_device_ptr( DEV_PTRS_K06 )
{
			for (int64_t inucl = 0; inucl < nucleus_num; ++inucl) {
				nucleus_range[inucl] = 0.;
				for (int64_t ishell = nucleus_index[inucl];
					 ishell <
					 nucleus_index[inucl] + nucleus_shell_num[inucl];
					 ++ishell) {
					for (int64_t iprim = shell_prim_index[ishell];
						 iprim < shell_prim_index[ishell] +
									 shell_prim_num[ishell];
						 ++iprim) {
						double range = 1. / exponent[iprim];
						nucleus_range[inucl] =
							nucleus_range[inucl] > range
								? nucleus_range[inucl]
								: range;
					}
				}
			}
}
		}
	}

	rc = qmckl_finalize_nucleus_basis_hpc_device(context);

	return rc;
}

/* Provided check  */

bool qmckl_nucleus_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->nucleus.provided;
}

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num) {
	int32_t mask = 1 << 0;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	if (num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_nucleus_num", "num <= 0");
	}

	ctx->nucleus.num = num;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_nucleus_num_device(const qmckl_context_device context,
							 int64_t *const num) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_nucleus_num",
									 "num is a null pointer");
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->nucleus.uninitialized & mask) != 0) {
		*num = (int64_t)0;
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_nucleus_num",
									 "nucleus data is not provided");
	}

	assert(ctx->nucleus.num >= (int64_t)0);
	*num = ctx->nucleus.num;

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_nucleus_coord_device(qmckl_context_device context, char transp,
							   double *coord, int64_t size_max) {
	int32_t mask = 1 << 2;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;

	if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
		return qmckl_failwith_device(context, QMCKL_ALREADY_SET_DEVICE,
									 "qmckl_set_nucleus_*", NULL);
	}

	qmckl_exit_code_device rc;

	const int64_t nucl_num = (int64_t)ctx->nucleus.num;

	if (ctx->nucleus.coord.data != NULL) {
		rc = qmckl_matrix_free_device(context, &(ctx->nucleus.coord));
		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	ctx->nucleus.coord = qmckl_matrix_alloc_device(context, nucl_num, 3);

	if (ctx->nucleus.coord.data == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_nucleus_coord", NULL);
	}

	if (size_max < 3 * nucl_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_4_DEVICE,
									 "qmckl_set_nucleus_coord",
									 "Array too small");
	}

	if (transp == 'N') {
		qmckl_matrix_device At;
		At = qmckl_matrix_alloc_device(context, 3, nucl_num);
		rc = qmckl_matrix_of_double_device(context, coord, 3 * nucl_num, &At);
		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
		rc = qmckl_transpose_device(context, At, ctx->nucleus.coord);
	} else {
		rc = qmckl_matrix_of_double_device(context, coord, nucl_num * 3,
										   &(ctx->nucleus.coord));
	}
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}

/* Sets the nuclear charges of all the atoms. */
qmckl_exit_code_device
qmckl_set_nucleus_charge_device(qmckl_context_device context, double *charge,
								int64_t size_max) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_NULL_CONTEXT_DEVICE,
									 "qmckl_set_nucleus_charge_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (charge == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "charge is a null pointer");
	}

	int64_t num;
	qmckl_exit_code_device rc;

	rc = qmckl_get_nucleus_num_device(context, &num);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	if (num > size_max) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "Array too small");
	}

	ctx->nucleus.charge = qmckl_vector_alloc_device(context, num);
	rc = qmckl_vector_of_double_device(context, charge, num,
									   &(ctx->nucleus.charge));

	if (rc != QMCKL_SUCCESS_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_nucleus_charge_device",
									 "Error in vector->double* conversion");
	}

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS_DEVICE;
}
