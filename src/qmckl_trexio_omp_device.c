#include "qmckl_trexio_device.h"

// This file provides wrappers to standard QMCkl functions accessible with the
// _device suffix. Includes OpenMP implementations.
// All getters/setters work with device pointers. Except when initializing the
// context from a TREXIO file, transferring data from/to the GPU must be handled
// by the user.

//**********
// ELECTRON/POINT GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_point_device(qmckl_context_device context,
									   const char transp, const int64_t num,
									   const double *coord,
									   const int64_t size_max) {

	const size_t device_id = qmckl_get_device_id(context);
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	if (size_max < 3 * num) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_4,
							  "qmckl_set_point_device", "Array too small");
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_point_device",
							  "transp should be 'N' or 'T'");
	}

	if (coord == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_point_device",
							  "coord is a NULL pointer");
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	qmckl_exit_code rc;
	if (num != ctx->point.num) {

		if (ctx->point.coord.data != NULL) {
			rc = qmckl_matrix_free_device(context, &(ctx->point.coord));
			assert(rc == QMCKL_SUCCESS);
		}

		ctx->point.coord = qmckl_matrix_alloc_device(context, num, 3);
		if (ctx->point.coord.data == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED, "qmckl_set_point",
								  NULL);
		}
	};

	ctx->point.num = num;

	double *a = ctx->point.coord.data;
	int size_0 = ctx->point.coord.size[0];
	if (transp == 'T') {
#pragma omp target is_device_ptr(a, coord)
		{
			for (int64_t i = 0; i < 3 * num; ++i) {
				a[i] = coord[i];
			}
		}
	} else {
#pragma omp target is_device_ptr(a, coord)
		{
			for (int64_t i = 0; i < num; ++i) {
				a[i] = coord[3 * i];
				a[i + size_0] = coord[3 * i + 1];
				a[i + 2 * size_0] = coord[3 * i + 1];
			}
		}
	}

	/* Increment the date of the context */
	rc = qmckl_context_touch_device(context);
	assert(rc == QMCKL_SUCCESS);

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_electron_coord_device(qmckl_context context,
												const char transp,
												const int64_t walk_num,
												const double *coord,
												const int64_t size_max) {

	const size_t device_id = qmckl_get_device_id(context);
	int32_t mask = 0; // coord can be changed

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
		return qmckl_failwith(context, QMCKL_ALREADY_SET,
							  "qmckl_set_electron_*", NULL);
	}

	if (transp != 'N' && transp != 'T') {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_electron_coord_device",
							  "transp should be 'N' or 'T'");
	}

	if (walk_num <= 0) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_electron_coord_device",
							  "walk_num <= 0");
	}

	if (coord == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_electron_coord_device",
							  "coord is a null pointer");
	}

	const int64_t elec_num = ctx->electron.num;

	if (elec_num == 0L) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_electron_coord_device",
							  "elec_num is not set");
	}

	/* Swap pointers */
	qmckl_walker tmp = ctx->electron.walker_old;
	ctx->electron.walker_old = ctx->electron.walker;
	ctx->electron.walker = tmp;

	memcpy(&(ctx->point), &(ctx->electron.walker.point),
		   sizeof(qmckl_point_struct));

	qmckl_exit_code rc;
	rc = qmckl_set_point_device(context, transp, walk_num * elec_num, coord,
								size_max);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->electron.walker.num = walk_num;
	memcpy(&(ctx->electron.walker.point), &(ctx->point),
		   sizeof(qmckl_point_struct));

	return QMCKL_SUCCESS;
}

//**********
// NUCLEUS GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_nucleus_charge_device(qmckl_context_device context,
												const double *charge,
												const int64_t size_max) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_charge_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (charge == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_nucleus_charge_device",
							  "charge is a null pointer");
	}

	int64_t num;
	qmckl_exit_code rc;

	rc = qmckl_get_nucleus_num(context, &num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	if (num > size_max) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_nucleus_charge_device",
							  "Array too small");
	}

	ctx->nucleus.charge = qmckl_vector_alloc_device(context, num);
	rc = qmckl_vector_of_double_device(context, charge, num,
									   &(ctx->nucleus.charge));

	if (rc != QMCKL_SUCCESS) {
		return qmckl_failwith(context, QMCKL_FAILURE,
							  "qmckl_set_nucleus_charge_device",
							  "Error in vector->double* conversion");
	}

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_nucleus_coord_device(qmckl_context_device context,
											   const char transp,
											   const double *coord,
											   const int64_t size_max) {
	int32_t mask = 1 << 2;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NULL_CONTEXT,
							  "qmckl_set_nucleus_coord_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	qmckl_exit_code rc;

	const int64_t nucl_num = (int64_t)ctx->nucleus.num;

	if (ctx->nucleus.coord.data != NULL) {
		rc = qmckl_matrix_free_device(context, &(ctx->nucleus.coord));
		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	ctx->nucleus.coord = qmckl_matrix_alloc_device(context, nucl_num, 3);

#pragma use_device_ptr(ctx->nucleus.coord.data)
	{
		if (ctx->nucleus.coord.data == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_set_nucleus_coord_device", NULL);
		}
	}

	if (size_max < 3 * nucl_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_4,
							  "qmckl_set_nucleus_coord_device",
							  "Array too small");
	}

	if (transp == 'N') {
		qmckl_matrix At;
		At = qmckl_matrix_alloc_device(context, 3, nucl_num);
		rc = qmckl_matrix_of_double_device(context, coord, 3 * nucl_num, &At);
		if (rc != QMCKL_SUCCESS)
			return rc;
		rc = qmckl_transpose_device(context, At, ctx->nucleus.coord);
	} else {
		rc = qmckl_matrix_of_double_device(context, coord, nucl_num * 3,
										   &(ctx->nucleus.coord));
	}
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->nucleus.uninitialized &= ~mask;
	ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context) {

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

	int device_id = qmckl_get_device_id(context);

	mem_info.size = ctx->nucleus.num * sizeof(int32_t);
	ctx->ao_basis.prim_num_per_nucleus =
		(int32_t *)qmckl_malloc_device(context, mem_info);

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

#pragma omp target map(tofrom                                                  \
					   : shell_max_ptr[:1], prim_max_ptr                       \
					   [:1])                                                   \
	is_device_ptr(nucleus_shell_num, nucleus_index, shell_prim_num,            \
				  prim_num_per_nucleus)
	{

		for (int inucl = 0; inucl < nucl_num; ++inucl) {
			shell_max_ptr[0] = nucleus_shell_num[inucl] > shell_max_ptr[0]
								   ? nucleus_shell_num[inucl]
								   : shell_max_ptr[0];

			int64_t prim_num = 0;
			const int64_t ishell_start = nucleus_index[inucl];
			const int64_t ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl];
			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {
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
	double *expo_expo = omp_target_alloc(prim_max * sizeof(double), device_id);
	int64_t *expo_index =
		omp_target_alloc(prim_max * sizeof(double), device_id);

	double *coef =
		omp_target_alloc(shell_max * prim_max * sizeof(double), device_id);
	double *newcoef = omp_target_alloc(prim_max * sizeof(double), device_id);

	int64_t *newidx = omp_target_alloc(prim_max * sizeof(int64_t), device_id);

	int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;
	double *exponent = ctx->ao_basis.exponent;
	double *coefficient_normalized = ctx->ao_basis.coefficient_normalized;

	double *expo_per_nucleus_data = ctx->ao_basis.expo_per_nucleus.data;
	int expo_per_nucleus_s0 = ctx->ao_basis.expo_per_nucleus.size[0];

	double *coef_per_nucleus_data = ctx->ao_basis.coef_per_nucleus.data;
	int coef_per_nucleus_s0 = ctx->ao_basis.coef_per_nucleus.size[0];
	int coef_per_nucleus_s1 = ctx->ao_basis.coef_per_nucleus.size[1];

#pragma omp target is_device_ptr(                                              \
	expo_expo, expo_index, coef, newcoef, nucleus_index, shell_prim_index,     \
	nucleus_shell_num, exponent, coefficient_normalized, shell_prim_num,       \
	expo_per_nucleus_data, coef_per_nucleus_data, prim_num_per_nucleus,        \
	newidx)
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

			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {

				const int64_t iprim_start = shell_prim_index[ishell];
				const int64_t iprim_end =
					shell_prim_index[ishell] + shell_prim_num[ishell];
				for (int64_t iprim = iprim_start; iprim < iprim_end; ++iprim) {
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
			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {

				for (int i = 0; i < prim_max; i++) {
					newcoef[i] = 0;
				}
				const int64_t iprim_start = shell_prim_index[ishell];
				const int64_t iprim_end =
					shell_prim_index[ishell] + shell_prim_num[ishell];

				for (int64_t iprim = iprim_start; iprim < iprim_end; ++iprim) {
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
					coef_per_nucleus_data[(i) + coef_per_nucleus_s0 *
													((j) + coef_per_nucleus_s1 *
															   (inucl))] =
						coef[j * prim_max + i];
				}
			}
		}
	}
	// End of target region

	omp_target_free(expo_expo, device_id);
	omp_target_free(expo_index, device_id);
	omp_target_free(coef, device_id);
	omp_target_free(newcoef, device_id);
	omp_target_free(newidx, device_id);

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_finalize_nucleus_basis_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_finalize_nucleus_basis_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t nucl_num = 0;

	qmckl_exit_code rc = qmckl_get_nucleus_num(context, &nucl_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* nucleus_prim_index */
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = (ctx->nucleus.num + (int64_t)1) * sizeof(int64_t);

		ctx->ao_basis.nucleus_prim_index =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.nucleus_prim_index == NULL) {
			return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
								  "ao_basis.nucleus_prim_index", NULL);
		}

		// Extract arrays from context
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_prim_index = ctx->ao_basis.nucleus_prim_index;
		int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;

		int prim_num = ctx->ao_basis.prim_num;

#pragma omp target is_device_ptr(nucleus_index, nucleus_prim_index,            \
								 shell_prim_index)
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
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ctx->ao_basis.prim_num * sizeof(double);

		ctx->ao_basis.coefficient_normalized =
			(double *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.coefficient_normalized == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
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

#pragma omp target is_device_ptr(shell_prim_index, shell_prim_num,             \
								 coefficient_normalized, coefficient,          \
								 prim_factor, shell_factor)
		{
			for (int64_t ishell = 0; ishell < shell_num; ++ishell) {
				for (int64_t iprim = shell_prim_index[ishell];
					 iprim < shell_prim_index[ishell] + shell_prim_num[ishell];
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
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ctx->nucleus.num * sizeof(int32_t);

		ctx->ao_basis.nucleus_max_ang_mom =
			(int32_t *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.nucleus_max_ang_mom == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "ao_basis.nucleus_max_ang_mom", NULL);
		}

		// Extract arrays from context
		int32_t *nucleus_max_ang_mom = ctx->ao_basis.nucleus_max_ang_mom;
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_shell_num = ctx->ao_basis.nucleus_shell_num;
		int32_t *shell_ang_mom = ctx->ao_basis.shell_ang_mom;

#pragma omp target is_device_ptr(nucleus_max_ang_mom, nucleus_index,           \
								 nucleus_shell_num, shell_ang_mom)
		{
#pragma omp parallel for
			for (int64_t inucl = 0; inucl < nucl_num; ++inucl) {
				nucleus_max_ang_mom[inucl] = 0;
				for (int64_t ishell = nucleus_index[inucl];
					 ishell < nucleus_index[inucl] + nucleus_shell_num[inucl];
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
			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size = ctx->nucleus.num * sizeof(double);

			ctx->ao_basis.nucleus_range =
				(double *)qmckl_malloc_device(context, mem_info);

			if (ctx->ao_basis.nucleus_range == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
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

#pragma omp target is_device_ptr(nucleus_range, nucleus_index,                 \
								 nucleus_shell_num, shell_prim_index,          \
								 shell_prim_num, exponent)
			{
				for (int64_t inucl = 0; inucl < nucleus_num; ++inucl) {
					nucleus_range[inucl] = 0.;
					for (int64_t ishell = nucleus_index[inucl];
						 ishell <
						 nucleus_index[inucl] + nucleus_shell_num[inucl];
						 ++ishell) {
						for (int64_t iprim = shell_prim_index[ishell];
							 iprim <
							 shell_prim_index[ishell] + shell_prim_num[ishell];
							 ++iprim) {
							double range = 1. / exponent[iprim];
							nucleus_range[inucl] = nucleus_range[inucl] > range
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

//**********
// AO GETTERS/SETTERS
//**********

qmckl_exit_code
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context) {

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

	int device_id = qmckl_get_device_id(context);

	mem_info.size = ctx->nucleus.num * sizeof(int32_t);
	ctx->ao_basis.prim_num_per_nucleus =
		(int32_t *)qmckl_malloc_device(context, mem_info);

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

#pragma omp target map(tofrom                                                  \
					   : shell_max_ptr[:1], prim_max_ptr                       \
					   [:1])                                                   \
	is_device_ptr(nucleus_shell_num, nucleus_index, shell_prim_num,            \
				  prim_num_per_nucleus)
	{

		for (int inucl = 0; inucl < nucl_num; ++inucl) {
			shell_max_ptr[0] = nucleus_shell_num[inucl] > shell_max_ptr[0]
								   ? nucleus_shell_num[inucl]
								   : shell_max_ptr[0];

			int64_t prim_num = 0;
			const int64_t ishell_start = nucleus_index[inucl];
			const int64_t ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl];
			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {
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
	double *expo_expo = omp_target_alloc(prim_max * sizeof(double), device_id);
	int64_t *expo_index =
		omp_target_alloc(prim_max * sizeof(double), device_id);

	double *coef =
		omp_target_alloc(shell_max * prim_max * sizeof(double), device_id);
	double *newcoef = omp_target_alloc(prim_max * sizeof(double), device_id);

	int64_t *newidx = omp_target_alloc(prim_max * sizeof(int64_t), device_id);

	int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;
	double *exponent = ctx->ao_basis.exponent;
	double *coefficient_normalized = ctx->ao_basis.coefficient_normalized;

	double *expo_per_nucleus_data = ctx->ao_basis.expo_per_nucleus.data;
	int expo_per_nucleus_s0 = ctx->ao_basis.expo_per_nucleus.size[0];

	double *coef_per_nucleus_data = ctx->ao_basis.coef_per_nucleus.data;
	int coef_per_nucleus_s0 = ctx->ao_basis.coef_per_nucleus.size[0];
	int coef_per_nucleus_s1 = ctx->ao_basis.coef_per_nucleus.size[1];

#pragma omp target is_device_ptr(                                              \
	expo_expo, expo_index, coef, newcoef, nucleus_index, shell_prim_index,     \
	nucleus_shell_num, exponent, coefficient_normalized, shell_prim_num,       \
	expo_per_nucleus_data, coef_per_nucleus_data, prim_num_per_nucleus,        \
	newidx)
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

			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {

				const int64_t iprim_start = shell_prim_index[ishell];
				const int64_t iprim_end =
					shell_prim_index[ishell] + shell_prim_num[ishell];
				for (int64_t iprim = iprim_start; iprim < iprim_end; ++iprim) {
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
			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {

				for (int i = 0; i < prim_max; i++) {
					newcoef[i] = 0;
				}
				const int64_t iprim_start = shell_prim_index[ishell];
				const int64_t iprim_end =
					shell_prim_index[ishell] + shell_prim_num[ishell];

				for (int64_t iprim = iprim_start; iprim < iprim_end; ++iprim) {
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
					coef_per_nucleus_data[(i) + coef_per_nucleus_s0 *
													((j) + coef_per_nucleus_s1 *
															   (inucl))] =
						coef[j * prim_max + i];
				}
			}
		}
	}
	// End of target region

	omp_target_free(expo_expo, device_id);
	omp_target_free(expo_index, device_id);
	omp_target_free(coef, device_id);
	omp_target_free(newcoef, device_id);
	omp_target_free(newidx, device_id);

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finalize_ao_basis_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_finalize_ao_basis_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int64_t nucl_num = 0;

	qmckl_exit_code rc = qmckl_get_nucleus_num(context, &nucl_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* nucleus_prim_index */
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = (ctx->nucleus.num + (int64_t)1) * sizeof(int64_t);

		ctx->ao_basis.nucleus_prim_index =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.nucleus_prim_index == NULL) {
			return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
								  "ao_basis.nucleus_prim_index", NULL);
		}

		// Extract arrays from context
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_prim_index = ctx->ao_basis.nucleus_prim_index;
		int64_t *shell_prim_index = ctx->ao_basis.shell_prim_index;

		int prim_num = ctx->ao_basis.prim_num;

#pragma omp target is_device_ptr(nucleus_index, nucleus_prim_index,            \
								 shell_prim_index)
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
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ctx->ao_basis.prim_num * sizeof(double);

		ctx->ao_basis.coefficient_normalized =
			(double *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.coefficient_normalized == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
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

#pragma omp target is_device_ptr(shell_prim_index, shell_prim_num,             \
								 coefficient_normalized, coefficient,          \
								 prim_factor, shell_factor)
		{
			for (int64_t ishell = 0; ishell < shell_num; ++ishell) {
				for (int64_t iprim = shell_prim_index[ishell];
					 iprim < shell_prim_index[ishell] + shell_prim_num[ishell];
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
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ctx->nucleus.num * sizeof(int32_t);

		ctx->ao_basis.nucleus_max_ang_mom =
			(int32_t *)qmckl_malloc_device(context, mem_info);

		if (ctx->ao_basis.nucleus_max_ang_mom == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "ao_basis.nucleus_max_ang_mom", NULL);
		}

		// Extract arrays from context
		int32_t *nucleus_max_ang_mom = ctx->ao_basis.nucleus_max_ang_mom;
		int64_t *nucleus_index = ctx->ao_basis.nucleus_index;
		int64_t *nucleus_shell_num = ctx->ao_basis.nucleus_shell_num;
		int32_t *shell_ang_mom = ctx->ao_basis.shell_ang_mom;

#pragma omp target is_device_ptr(nucleus_max_ang_mom, nucleus_index,           \
								 nucleus_shell_num, shell_ang_mom)
		{
#pragma omp parallel for
			for (int64_t inucl = 0; inucl < nucl_num; ++inucl) {
				nucleus_max_ang_mom[inucl] = 0;
				for (int64_t ishell = nucleus_index[inucl];
					 ishell < nucleus_index[inucl] + nucleus_shell_num[inucl];
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
			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size = ctx->nucleus.num * sizeof(double);

			ctx->ao_basis.nucleus_range =
				(double *)qmckl_malloc_device(context, mem_info);

			if (ctx->ao_basis.nucleus_range == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
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

#pragma omp target is_device_ptr(nucleus_range, nucleus_index,                 \
								 nucleus_shell_num, shell_prim_index,          \
								 shell_prim_num, exponent)
			{
				for (int64_t inucl = 0; inucl < nucleus_num; ++inucl) {
					nucleus_range[inucl] = 0.;
					for (int64_t ishell = nucleus_index[inucl];
						 ishell <
						 nucleus_index[inucl] + nucleus_shell_num[inucl];
						 ++ishell) {
						for (int64_t iprim = shell_prim_index[ishell];
							 iprim <
							 shell_prim_index[ishell] + shell_prim_num[ishell];
							 ++iprim) {
							double range = 1. / exponent[iprim];
							nucleus_range[inucl] = nucleus_range[inucl] > range
													   ? nucleus_range[inucl]
													   : range;
						}
					}
				}
			}
		}
	}

	rc = qmckl_finalize_ao_basis_hpc_device(context);

	return rc;
}

qmckl_exit_code qmckl_set_ao_basis_type_device(qmckl_context_device context,
											   const char basis_type) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_type_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (basis_type != 'G' && basis_type != 'S') {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_type_device", NULL);
	}

	int32_t mask = 1;

	ctx->ao_basis.type = basis_type;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									const int64_t shell_num) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_num_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (shell_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num <= 0");
	}

	const int64_t prim_num = ctx->ao_basis.prim_num;

	if (0L < prim_num && prim_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num > prim_num");
	}

	int32_t mask = 1 << 1;

	ctx->ao_basis.shell_num = shell_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
												   const int64_t prim_num) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_num_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (prim_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_prim_num_device",
							  "prim_num must be positive");
	}

	const int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_prim_num_device",
							  "shell_num is not set");
	}

	if (prim_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_prim_num_device",
							  "prim_num < shell_num");
	}

	int32_t mask = 1 << 2;

	ctx->ao_basis.prim_num = prim_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_ao_num_device(qmckl_context_device context,
												 const int64_t ao_num) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_ao_num_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (ao_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "ao_num must be positive");
	}

	const int64_t shell_num = ctx->ao_basis.shell_num;
	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "shell_num is not set");
	}

	if (ao_num < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_ao_basis_shell_num_device",
							  "ao_num < shell_num");
	}

	int32_t mask = 1 << 12;

	ctx->ao_basis.ao_num = ao_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_index_device(qmckl_context_device context,
										const int64_t *nucleus_index,
										const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_nucleus_index_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 4;

	const int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_index_device",
							  "nucl_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_index_device",
							  "input array too small");
	}

	if (ctx->ao_basis.nucleus_index != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_index);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_nucleus_index_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = nucl_num * sizeof(int64_t);
	int64_t *new_array = (int64_t *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_nucleus_index_device", NULL);
	}

	omp_target_memcpy(new_array, nucleus_index, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.nucleus_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											const int64_t *nucleus_shell_num,
											const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 3;

	const int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  "shell_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  "input array too small");
	}

	if (ctx->ao_basis.nucleus_shell_num != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_shell_num);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith(context, rc,
								  "qmckl_set_ao_basis_nucleus_shell_num_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = nucl_num * sizeof(int64_t);
	int64_t *new_array = (int64_t *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_nucleus_shell_num_device",
							  NULL);
	}

	omp_target_memcpy(new_array, nucleus_shell_num, mem_info.size, 0, 0,
					  device_id, device_id);

	ctx->ao_basis.nucleus_shell_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_ang_mom_device(qmckl_context_device context,
										const int32_t *shell_ang_mom,
										const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 5;

	const int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num == 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_ang_mom_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_ang_mom_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_ang_mom != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_ang_mom);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_ang_mom_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = shell_num * sizeof(int32_t);
	int32_t *new_array = (int32_t *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
	}

	omp_target_memcpy(new_array, shell_ang_mom, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.shell_ang_mom = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_num_device(qmckl_context_device context,
										 const int64_t *shell_prim_num,
										 const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(
			(qmckl_context)context, QMCKL_INVALID_CONTEXT,
			"qmckl_set_ao_basis_nucleus_shell_prim_num_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 6;

	const int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_prim_num_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_prim_num_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_prim_num != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_num);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_prim_num_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = shell_num * sizeof(int64_t);
	int64_t *new_array = (int64_t *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_prim_num_device", NULL);
	}

	omp_target_memcpy(new_array, shell_prim_num, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.shell_prim_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_index_device(qmckl_context_device context,
										   const int64_t *shell_prim_index,
										   const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 7;

	const int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_prim_index != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_index);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_prim_index_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = shell_num * sizeof(int64_t);
	int64_t *new_array = (int64_t *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_prim_index_device",
							  NULL);
	}

	omp_target_memcpy(new_array, shell_prim_index, mem_info.size, 0, 0,
					  device_id, device_id);

	ctx->ao_basis.shell_prim_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   const double *shell_factor,
									   const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_shell_factor_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 8;

	const int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_shell_factor_device",
							  "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_shell_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.shell_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.shell_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_shell_factor_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = shell_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_shell_factor_device", NULL);
	}

	omp_target_memcpy(new_array, shell_factor, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.shell_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
												   const double *exponent,
												   const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_exponent_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 9;

	const int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_exponent_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_exponent_device",
							  "input array too small");
	}

	if (ctx->ao_basis.exponent != NULL) {
		qmckl_exit_code rc = qmckl_free_device(context, ctx->ao_basis.exponent);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_exponent_device", NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = prim_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_exponent_device", NULL);
	}

	omp_target_memcpy(new_array, exponent, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.exponent = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  const double *coefficient,
									  const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_coefficient_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 10;

	const int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_coefficient_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_coefficient_device",
							  "input array too small");
	}

	if (ctx->ao_basis.coefficient != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.coefficient);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_coefficient_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = prim_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_coefficient_device", NULL);
	}

	omp_target_memcpy(new_array, coefficient, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.coefficient = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_ao_basis_prim_factor_device(qmckl_context context,
													  const double *prim_factor,
													  const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	// Accepts an host array and copies it on device

	int32_t mask = 1 << 11;

	const int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_prim_factor_device",
							  "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_prim_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.prim_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.prim_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_prim_factor_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = prim_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}

	omp_target_memcpy(new_array, prim_factor, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.prim_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									const double *ao_factor,
									const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_set_ao_basis_prim_factor_device", NULL);
	}
	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 13;

	const int64_t ao_num = ctx->ao_basis.ao_num;

	if (ao_num <= 0L) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_set_ao_basis_ao_factor_device",
							  "ao_num is not set");
	}

	if (size_max < ao_num) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_set_ao_basis_ao_factor_device",
							  "input array too small");
	}

	if (ctx->ao_basis.ao_factor != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->ao_basis.ao_factor);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_ao_basis_ao_factor_device", NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = ao_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);

	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_ao_basis_ao_factor_device", NULL);
	}

	omp_target_memcpy(new_array, ao_factor, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->ao_basis.ao_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

//**********
// MO GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_finalize_mo_basis_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_finalize_mo_basis_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size =
		ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);
	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_finalize_mo_basis_device", NULL);
	}

	assert(ctx->mo_basis.coefficient != NULL);

	if (ctx->mo_basis.coefficient_t != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->mo_basis.coefficient_t);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_finalize_mo_basis_device", NULL);
		}
	}

	double *coefficient = ctx->mo_basis.coefficient;

	int64_t ao_num = ctx->ao_basis.ao_num;
	int64_t mo_num = ctx->mo_basis.mo_num;

#pragma omp target is_device_ptr(new_array, coefficient)
	{
#pragma omp parallel for collapse(2)
		for (int64_t i = 0; i < ao_num; ++i) {
			for (int64_t j = 0; j < mo_num; ++j) {
				new_array[i * mo_num + j] = coefficient[j * ao_num + i];
			}
		}
	}

	ctx->mo_basis.coefficient_t = new_array;
	qmckl_exit_code rc = QMCKL_SUCCESS;
	return rc;
}

qmckl_exit_code qmckl_set_mo_basis_mo_num_device(qmckl_context_device context,
												 const int64_t mo_num) {

	int32_t mask = 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALREADY_SET,
							  "qmckl_set_mo_basis_mo_num_device", NULL);
	}

	if (mo_num <= 0) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_set_mo_basis_mo_num_device",
							  "mo_num <= 0");
	}

	ctx->mo_basis.mo_num = mo_num;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_mo_basis_coefficient_device(qmckl_context context,
									  const double *coefficient) {

	int32_t mask = 1 << 1;

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int device_id = qmckl_get_device_id(context);

	if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALREADY_SET,
							  "qmckl_set_mo_basis_coefficient_device", NULL);
	}

	if (ctx->mo_basis.coefficient != NULL) {
		qmckl_exit_code rc =
			qmckl_free_device(context, ctx->mo_basis.coefficient);
		if (rc != QMCKL_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, rc,
								  "qmckl_set_mo_basis_coefficient_device",
								  NULL);
		}
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size =
		ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
	double *new_array = (double *)qmckl_malloc_device(context, mem_info);
	if (new_array == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_set_mo_basis_coefficient_device", NULL);
	}

	omp_target_memcpy(new_array, coefficient, mem_info.size, 0, 0, device_id,
					  device_id);

	ctx->mo_basis.coefficient = new_array;

	ctx->mo_basis.uninitialized &= ~mask;
	ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
	if (ctx->mo_basis.provided) {
		qmckl_exit_code rc_ = qmckl_finalize_mo_basis_device(context);
		if (rc_ != QMCKL_SUCCESS)
			return rc_;
	}

	return QMCKL_SUCCESS;
}

//**********
// CONTEXT FILL
//**********

qmckl_exit_code
qmckl_trexio_read_electron_X_device(qmckl_context_device context,
									trexio_t *const file) {

	printf("[qmckl_trexio_read_electron_X_device] In\n");
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	int rcio = 0;

	int64_t up_num = 0L;
	int64_t dn_num = 0L;

	rcio = trexio_read_electron_up_num_64(file, &up_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_electron_up_num",
							  trexio_string_of_error(rcio));
	}

	assert(up_num >= 0L);

	rcio = trexio_read_electron_dn_num_64(file, &dn_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_electron_dn_num",
							  trexio_string_of_error(rcio));
	}

	assert(dn_num >= 0L);

	qmckl_exit_code rc;
	rc = qmckl_set_electron_num_device(context, up_num, dn_num);
	printf("[qmckl_trexio_read_electron_X_device] Out\n");
	return rc;
}

qmckl_exit_code qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
												   trexio_t *const file) {
	printf("[qmckl_trexio_read_nucleus_X_device] In\n");
	assert(context != (qmckl_context)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;

	int64_t nucleus_num = 0L;

	rcio = trexio_read_nucleus_num_64(file, &nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_nucleus_num",
							  trexio_string_of_error(rcio));
	}

	assert(nucleus_num > 0);
	rc = qmckl_set_nucleus_num_device(context, nucleus_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = nucleus_num * sizeof(double);

		double *nucl_charge_h = (double *)qmckl_malloc_host(context, mem_info);
		double *nucl_charge_d =
			(double *)qmckl_malloc_device(context, mem_info);

		if (nucl_charge_h == NULL || nucl_charge_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_nucleus_X_device", NULL);
		}

		assert(nucl_charge_h != NULL && nucl_charge_d != NULL);

		rcio = trexio_read_safe_nucleus_charge_64(file, nucl_charge_h,
												  nucleus_num);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_nucleus_charge",
								  trexio_string_of_error(rcio));
		}

		qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge_h, mem_info.size);
		rc = qmckl_set_nucleus_charge_device(context, nucl_charge_d,
											 nucleus_num);

		qmckl_free_host(context, nucl_charge_h);
		qmckl_free_device(context, nucl_charge_d);
		nucl_charge_h = NULL;
		nucl_charge_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
	mem_info.size = nucleus_num * 3 * sizeof(double);

	double *nucl_coord_h = (double *)qmckl_malloc_host(context, mem_info);
	double *nucl_coord_d = (double *)qmckl_malloc_device(context, mem_info);

	if (nucl_coord_h == NULL || nucl_coord_d == NULL) {
		return qmckl_failwith((qmckl_context)context, QMCKL_ALLOCATION_FAILED,
							  "qmckl_trexio_read_nucleus_X_device", NULL);
	}

	assert(nucl_coord_h != NULL && nucl_coord_d != NULL);

	rcio =
		trexio_read_safe_nucleus_coord_64(file, nucl_coord_h, 3 * nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_nucleus_charge",
							  trexio_string_of_error(rcio));
	}

	qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord_h, mem_info.size);
	rc = qmckl_set_nucleus_coord_device((qmckl_context)context, 'N',
										nucl_coord_d, 3 * nucleus_num);

	qmckl_free_host(context, nucl_coord_h);
	qmckl_free_device(context, nucl_coord_d);
	nucl_coord_h = NULL;
	nucl_coord_d = NULL;

	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	printf("[qmckl_trexio_read_nucleus_X_device] Out\n");
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_ao_X_device(qmckl_context context,
											  trexio_t *const file) {
	printf("[qmckl_trexio_read_ao_X_device] In\n");
	assert(context != (qmckl_context)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;
	int64_t nucleus_num = 0L;

	rc = qmckl_get_nucleus_num_device(context, &nucleus_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

#define MAX_STR_LEN 1024
	char basis_type[MAX_STR_LEN];

	printf("[qmckl_trexio_read_ao_X_device] 0\n");
	rcio = trexio_read_basis_type(file, basis_type, MAX_STR_LEN);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_type",
							  trexio_string_of_error(rcio));
	}

	if (basis_type[0] == 'G') {
		rc = qmckl_set_ao_basis_type_device(context, basis_type[0]);
	} else {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_type", "Invalid basis type");
	}

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t shell_num = 0L;

	printf("[qmckl_trexio_read_ao_X_device] 1\n");
	rcio = trexio_read_basis_shell_num_64(file, &shell_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_shell_num",
							  trexio_string_of_error(rcio));
	}

	assert(shell_num > 0);
	rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t prim_num = 0L;

	rcio = trexio_read_basis_prim_num_64(file, &prim_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_basis_prim_num",
							  trexio_string_of_error(rcio));
	}

	assert(prim_num > 0);
	rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t ao_num = 0LL;

	printf("[qmckl_trexio_read_ao_X_device] 2\n");
	rcio = trexio_read_ao_num_64(file, &ao_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_ao_num",
							  trexio_string_of_error(rcio));
	}

	assert(ao_num > 0);
	rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	printf("[qmckl_trexio_read_ao_X_device] 3\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;
		int64_t *nucleus_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_index_d =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (nucleus_index_h == NULL || nucleus_index_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "index_X_device",
								  NULL);
		}

		assert(nucleus_index_h != NULL && nucleus_index_d != NULL);

		/* Allocate temporary array */
		mem_info.size = shell_num * sizeof(int64_t);
		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "index_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read in the temporary array */
		rcio =
			trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_nucleus_index",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
													 nucleus_num);
		if (rc != QMCKL_SUCCESS) {
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return rc;
		}

		for (int i = shell_num - 1; i >= 0; --i) {
			const int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_index_h);
				qmckl_free_device(context, nucleus_index_d);
				nucleus_index_h = NULL;
				nucleus_index_d = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_nucleus_index",
									  "Irrelevant data in TREXIO file");
			}
			nucleus_index_h[k] = i;
		}

		qmckl_memcpy_H2D(context, nucleus_index_d, nucleus_index_h,
						 size_backup);

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
													 shell_num);

		qmckl_free_host(context, nucleus_index_h);
		qmckl_free_device(context, nucleus_index_d);
		nucleus_index_h = NULL;
		nucleus_index_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 4\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *nucleus_shell_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_shell_num_d =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (nucleus_shell_num_h == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "shell_num_X_device",
								  NULL);
		}

		assert(nucleus_shell_num_h != NULL);

		/* Allocate temporary array */
		mem_info.size = shell_num * sizeof(int64_t);
		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, nucleus_shell_num_h);
			qmckl_free_device(context, nucleus_shell_num_d);
			nucleus_shell_num_h = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_nucleus_"
								  "shell_num_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read in the temporary array */
		rcio =
			trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			qmckl_free_host(context, nucleus_shell_num_h);
			qmckl_free_device(context, nucleus_shell_num_d);
			nucleus_shell_num_h = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_nucleus_shell_num",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < nucleus_num; ++i) {
			nucleus_shell_num_h[i] = 0;
		}

		for (int i = 0; i < shell_num; ++i) {
			const int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_shell_num_h);
				qmckl_free_device(context, nucleus_shell_num_d);
				nucleus_shell_num_h = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_nucleus_shell_num",
									  "Irrelevant data in TREXIO file");
			}
			nucleus_shell_num_h[k] += 1;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, nucleus_shell_num_d, nucleus_shell_num_h,
						 size_backup);
		rc = qmckl_set_ao_basis_nucleus_shell_num_device(
			context, nucleus_shell_num_d, shell_num);

		qmckl_free_host(context, nucleus_shell_num_h);
		qmckl_free_device(context, nucleus_shell_num_d);
		nucleus_shell_num_h = NULL;
		nucleus_shell_num_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 5\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int32_t);

		int32_t *shell_ang_mom_h =
			(int32_t *)qmckl_malloc_host(context, mem_info);
		int32_t *shell_ang_mom_d =
			(int32_t *)qmckl_malloc_device(context, mem_info);

		if (shell_ang_mom_h == NULL || shell_ang_mom_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "ang_mom_X_device",
								  NULL);
		}

		assert(shell_ang_mom_h != NULL && shell_ang_mom_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_ang_mom_32(file, shell_ang_mom_h,
													   shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_ang_mom_h);
			qmckl_free_device(context, shell_ang_mom_d);
			shell_ang_mom_h = NULL;
			shell_ang_mom_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_ang_mom",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, shell_ang_mom_d, shell_ang_mom_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom_d,
													 shell_num);

		qmckl_free_host(context, shell_ang_mom_h);
		qmckl_free_device(context, shell_ang_mom_d);
		shell_ang_mom_h = NULL;
		shell_ang_mom_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 6\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_num_d =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (shell_prim_num_h == NULL || shell_prim_num_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_num_X_device",
								  NULL);
		}

		assert(shell_prim_num_h != NULL && shell_prim_num_d != NULL);

		/* Allocate temporary array */
		mem_info.size = prim_num * sizeof(int64_t);

		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			qmckl_free_host(context, shell_prim_num_h);
			qmckl_free_device(context, shell_prim_num_d);
			shell_prim_num_h = NULL;
			shell_prim_num_d = NULL;
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_num_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_prim_num_h);
			qmckl_free_device(context, shell_prim_num_d);
			shell_prim_num_h = NULL;
			shell_prim_num_d = NULL;
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_prim_num",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < shell_num; ++i) {
			shell_prim_num_h[i] = 0;
		}

		for (int i = 0; i < prim_num; ++i) {
			const int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				qmckl_free_host(context, shell_prim_num_h);
				qmckl_free_device(context, shell_prim_num_d);
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_shell_prim_num",
									  "Irrelevant data in TREXIO file");
			}
			shell_prim_num_h[k] += 1;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, shell_prim_num_d, shell_prim_num_h,
						 size_backup);
		rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num_d,
													  shell_num);

		qmckl_free_host(context, shell_prim_num_h);
		qmckl_free_device(context, shell_prim_num_d);
		shell_prim_num_h = NULL;
		shell_prim_num_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 7\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_index_d =
			(int64_t *)qmckl_malloc_device(context, mem_info);

		if (shell_prim_index_h == NULL || shell_prim_index_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_index_X_device",
								  NULL);
		}

		assert(shell_prim_index_h != NULL && shell_prim_index_d != NULL);

		/* Allocate temporary array */
		mem_info.size = prim_num * sizeof(int64_t);

		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_shell_"
								  "prim_index_X_device",
								  NULL);
		}

		assert(tmp_array != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_prim_index_h);
			qmckl_free_device(context, shell_prim_index_d);
			shell_prim_index_h = NULL;
			shell_prim_index_d = NULL;
			qmckl_free_host(context, tmp_array);
			tmp_array = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_prim_index",
								  trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = prim_num - 1; i >= 0; --i) {
			const int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, shell_prim_index_h);
				qmckl_free_device(context, shell_prim_index_d);
				shell_prim_index_h = NULL;
				shell_prim_index_d = NULL;
				return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
									  "trexio_read_basis_shell_prim_index",
									  "Irrelevant data in TREXIO file");
			}
			shell_prim_index_h[k] = i;
		}

		qmckl_free_host(context, tmp_array);
		tmp_array = NULL;

		/* Store data */
		qmckl_memcpy_H2D(context, shell_prim_index_d, shell_prim_index_h,
						 size_backup);
		rc = qmckl_set_ao_basis_shell_prim_index_device(
			context, shell_prim_index_d, shell_num);

		qmckl_free_host(context, shell_prim_index_h);
		qmckl_free_device(context, shell_prim_index_d);
		shell_prim_index_h = NULL;
		shell_prim_index_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 8\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(double);

		double *shell_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *shell_factor_d =
			(double *)qmckl_malloc_device(context, mem_info);

		if (shell_factor_h == NULL || shell_factor_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_shell_factor_X_device", NULL);
		}

		assert(shell_factor_h != NULL && shell_factor_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_shell_factor_64(file, shell_factor_h,
													  shell_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, shell_factor_h);
			qmckl_free_device(context, shell_factor_d);
			shell_factor_h = NULL;
			shell_factor_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_shell_factor",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, shell_factor_d, shell_factor_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor_d,
													shell_num);

		qmckl_free_host(context, shell_factor_h);
		qmckl_free_device(context, shell_factor_d);
		shell_factor_h = NULL;
		shell_factor_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 9\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *exponent_h = (double *)qmckl_malloc_host(context, mem_info);
		double *exponent_d = (double *)qmckl_malloc_device(context, mem_info);

		if (exponent_h == NULL || exponent_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_basis_exponent_X", NULL);
		}

		assert(exponent_h != NULL && exponent_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_exponent_64(file, exponent_h, prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, exponent_h);
			qmckl_free_device(context, exponent_d);
			exponent_h = NULL;
			exponent_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_exponent",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, exponent_d, exponent_h, mem_info.size);
		rc = qmckl_set_ao_basis_exponent_device(context, exponent_d, prim_num);

		qmckl_free_host(context, exponent_h);
		qmckl_free_device(context, exponent_d);
		exponent_h = NULL;
		exponent_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 10\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *coefficient_h = (double *)qmckl_malloc_host(context, mem_info);
		double *coefficient_d =
			(double *)qmckl_malloc_device(context, mem_info);

		if (coefficient_h == NULL || coefficient_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_coefficient_X_device", NULL);
		}

		assert(coefficient_h != NULL && coefficient_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_coefficient_64(file, coefficient_h,
													 prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, coefficient_h);
			qmckl_free_device(context, coefficient_d);
			coefficient_h = NULL;
			coefficient_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_coefficient",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, coefficient_d, coefficient_h, mem_info.size);
		rc = qmckl_set_ao_basis_coefficient_device(context, coefficient_d,
												   prim_num);

		qmckl_free_host(context, coefficient_h);
		qmckl_free_device(context, coefficient_d);
		coefficient_h = NULL;
		coefficient_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 11\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *prim_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *prim_factor_d =
			(double *)qmckl_malloc_device(context, mem_info);

		if (prim_factor_h == NULL || prim_factor_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_basis_prim_factor_X_device", NULL);
		}

		assert(prim_factor_h != NULL && prim_factor_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_basis_prim_factor_64(file, prim_factor_h,
													 prim_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, prim_factor_h);
			qmckl_free_device(context, prim_factor_d);
			prim_factor_h = NULL;
			prim_factor_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_basis_prim_factor",
								  trexio_string_of_error(rcio));
		}

		/* Read data */
		qmckl_memcpy_H2D(context, prim_factor_d, prim_factor_h, mem_info.size);
		rc = qmckl_set_ao_basis_prim_factor_device(context, prim_factor_d,
												   prim_num);

		qmckl_free_host(context, prim_factor_h);
		qmckl_free_device(context, prim_factor_d);
		prim_factor_h = NULL;
		prim_factor_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_read_ao_X_device] 12\n");
	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

		/* Allocate array for data */
		mem_info.size = ao_num * sizeof(double);

		double *ao_normalization_h =
			(double *)qmckl_malloc_host(context, mem_info);
		double *ao_normalization_d =
			(double *)qmckl_malloc_device(context, mem_info);

		if (ao_normalization_h == NULL || ao_normalization_d == NULL) {
			return qmckl_failwith(
				(qmckl_context)context, QMCKL_ALLOCATION_FAILED,
				"qmckl_trexio_read_ao_normalization_X_device", NULL);
		}

		assert(ao_normalization_h != NULL && ao_normalization_d != NULL);

		/* Read data */
		rcio = trexio_read_safe_ao_normalization_64(file, ao_normalization_h,
													ao_num);
		if (rcio != TREXIO_SUCCESS) {
			qmckl_free_host(context, ao_normalization_h);
			qmckl_free_device(context, ao_normalization_d);
			ao_normalization_h = NULL;
			ao_normalization_d = NULL;
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_ao_normalization",
								  trexio_string_of_error(rcio));
		}

		/* Store data */
		qmckl_memcpy_H2D(context, ao_normalization_d, ao_normalization_h,
						 mem_info.size);
		rc = qmckl_set_ao_basis_ao_factor_device(context, ao_normalization_d,
												 ao_num);

		qmckl_free_host(context, ao_normalization_h);
		qmckl_free_device(context, ao_normalization_d);
		ao_normalization_h = NULL;
		ao_normalization_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_ao_X_device] Out\n");
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_mo_X_device(qmckl_context_device context,
											  trexio_t *const file) {
	printf("[qmckl_trexio_mo_X_device] In\n");
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	qmckl_exit_code rc;
	int rcio = 0;
	int64_t ao_num = 0L;

	rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num);
	if (rc != QMCKL_SUCCESS)
		return rc;

	int64_t mo_num = 0L;

	rcio = trexio_read_mo_num_64(file, &mo_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "trexio_read_mo_num",
							  trexio_string_of_error(rcio));
	}

	assert(mo_num > 0);
	rc = qmckl_set_mo_basis_mo_num_device(context, mo_num);

	if (rc != QMCKL_SUCCESS)
		return rc;

	{
		qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
		mem_info.size = ao_num * mo_num * sizeof(double);

		double *mo_coef_h = (double *)qmckl_malloc_host(context, mem_info);
		double *mo_coef_d = (double *)qmckl_malloc_device(context, mem_info);

		if (mo_coef_h == NULL || mo_coef_d == NULL) {
			return qmckl_failwith((qmckl_context)context,
								  QMCKL_ALLOCATION_FAILED,
								  "qmckl_trexio_read_mo_X_device", NULL);
		}

		assert(mo_coef_h != NULL && mo_coef_d != NULL);

		rcio = trexio_read_mo_coefficient_64(file, mo_coef_h);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
								  "trexio_read_mo_coefficient",
								  trexio_string_of_error(rcio));
		}

		qmckl_memcpy_H2D(context, mo_coef_d, mo_coef_h, mem_info.size);
		rc = qmckl_set_mo_basis_coefficient_device(context, mo_coef_d);

		qmckl_free_host(context, mo_coef_h);
		qmckl_free_device(context, mo_coef_d);
		mo_coef_h = NULL;
		mo_coef_d = NULL;

		if (rc != QMCKL_SUCCESS)
			return rc;
	}

	printf("[qmckl_trexio_mo_X_device] Out\n");
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_trexio_read_device(const qmckl_context_device context,
										 const char *file_name,
										 const int64_t size_max) {
	printf("[qmckl_trexio_read_device] In\n");
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return false;
	}

	qmckl_exit_code rc;
	char file_name_new[size_max + 1];
	strncpy(file_name_new, file_name, size_max + 1);
	file_name_new[size_max] = '\0';

	trexio_t *file = qmckl_trexio_open_X(file_name_new, &rc);
	if (file == NULL) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_2,
							  "qmckl_trexio_read_device",
							  trexio_string_of_error(rc));
	}

	assert(file != NULL);

	rc = qmckl_trexio_read_electron_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device",
							  "Error reading electron");
	}

	rc = qmckl_trexio_read_nucleus_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device",
							  "Error reading nucleus");
	}

	rc = qmckl_trexio_read_ao_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_read_device", "Error reading AOs");
	}

	rc = qmckl_trexio_read_mo_X_device(context, file);
	if (rc != QMCKL_SUCCESS) {
		trexio_close(file);
		return qmckl_failwith((qmckl_context)context, QMCKL_FAILURE,
							  "qmckl_trexio_omp_read", "Error reading MOs");
	}

	trexio_close(file);
	file = NULL;
	printf("[qmckl_trexio_read_device] Out\n");
	return rc;
}
