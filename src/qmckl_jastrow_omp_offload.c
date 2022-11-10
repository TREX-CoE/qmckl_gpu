#include "include/qmckl_jastrow_gpu.h"

//**********
// COMPUTE
//**********

qmckl_exit_code qmckl_compute_tmp_c_omp_offload(
	const qmckl_context context, const int64_t cord_num, const int64_t elec_num,
	const int64_t nucl_num, const int64_t walk_num,
	const double *een_rescaled_e, const double *een_rescaled_n,
	double *const tmp_c) {

	if (context == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (cord_num <= 0) {
		return QMCKL_INVALID_ARG_2;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3;
	}

	if (nucl_num <= 0) {
		return QMCKL_INVALID_ARG_4;
	}

	// Compute array access strides:
	// For tmp_c...
	const int64_t stride_k_c = elec_num;
	const int64_t stride_j_c = stride_k_c * nucl_num;
	const int64_t stride_i_c = stride_j_c * (cord_num + 1);
	const int64_t stride_nw_c = stride_i_c * cord_num;
	// For een_rescaled_e...
	const int64_t stride_m_e = elec_num;
	const int64_t stride_i_e = stride_m_e * elec_num;
	const int64_t stride_nw_e = stride_i_e * (cord_num + 1);
	// For een_rescaled_n...
	const int64_t stride_k_n = elec_num;
	const int64_t stride_j_n = stride_k_n * nucl_num;
	const int64_t stride_nw_n = stride_j_n * (cord_num + 1);

	const int64_t size_tmp_c =
		elec_num * nucl_num * (cord_num + 1) * cord_num * walk_num;
	const int64_t size_e = walk_num * (cord_num + 1) * elec_num * elec_num;
	const int64_t size_n = walk_num * (cord_num + 1) * nucl_num * elec_num;

#pragma omp target data map(from                                               \
							: tmp_c [0:size_tmp_c])                            \
	map(to                                                                     \
		: een_rescaled_e [0:size_e], een_rescaled_n [0:size_n])
	{
#pragma omp target teams distribute parallel for collapse(4)
		for (int nw = 0; nw < walk_num; ++nw) {
			for (int i = 0; i < cord_num; ++i) {

				// Replacement for single DGEMM
				for (int jk = 0; jk < nucl_num * (cord_num + 1); jk++) {
					for (int l = 0; l < elec_num; l++) {

						int index_e_base =
							l + i * stride_i_e + nw * stride_nw_e;
						int index_n_base = jk * stride_k_n + nw * stride_nw_n;

						// Single reduction
						double sum = 0.;
#pragma omp simd reduction(+ : sum)
						for (int m = 0; m < elec_num; m++) {
							sum +=
								een_rescaled_e[index_e_base + m * stride_m_e] *
								een_rescaled_n[index_n_base + m];
						}
						tmp_c[l + jk * stride_k_c + i * stride_i_c +
							  nw * stride_nw_c] = sum;
					}
				}
			}
		}
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_dtmp_c_omp_offload(
	const qmckl_context context, const int64_t cord_num, const int64_t elec_num,
	const int64_t nucl_num, const int64_t walk_num,
	const double *een_rescaled_e_deriv_e, const double *een_rescaled_n,
	double *const dtmp_c) {

	if (context == QMCKL_NULL_CONTEXT) {
		return QMCKL_INVALID_CONTEXT;
	}

	if (cord_num <= 0) {
		return QMCKL_INVALID_ARG_2;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3;
	}

	if (nucl_num <= 0) {
		return QMCKL_INVALID_ARG_4;
	}

	// Compute strides...
	// For dtmp_c
	const int64_t stride_l_d = elec_num;
	const int64_t stride_k_d = stride_l_d * 4;
	const int64_t stride_j_d = stride_k_d * nucl_num;
	const int64_t stride_i_d = stride_j_d * (cord_num + 1);
	const int64_t stride_nw_d = stride_i_d * cord_num;
	// For een_rescaled_e_deriv_e
	const int64_t stride_l_e = elec_num;
	const int64_t stride_n_e = stride_l_e * 4;
	const int64_t stride_i_e = stride_n_e * elec_num;
	const int64_t stride_nw_e = stride_i_e * cord_num;
	// For een_rescaled_n
	const int64_t stride_k_n = elec_num;
	const int64_t stride_j_n = stride_k_n * nucl_num;
	const int64_t stride_nw_n = stride_j_n * (cord_num + 1);

	const int64_t size_dtmp_c =
		walk_num * cord_num * (cord_num + 1) * nucl_num * 4 * elec_num;
	const int64_t size_n = walk_num * (cord_num + 1) * nucl_num * elec_num;
	const int64_t size_e = walk_num * (cord_num + 1) * elec_num * 4 * elec_num;

	double sum = 0.;
#pragma omp target data map(from                                               \
							: dtmp_c [0:size_dtmp_c])                          \
	map(to                                                                     \
		: een_rescaled_e_deriv_e [0:size_e], een_rescaled_n [0:size_n])
	{
#pragma omp target teams distribute parallel for collapse(4)
		for (int nw = 0; nw < walk_num; nw++) {
			for (int i = 0; i < cord_num; i++) {

				// Single DGEMM
				for (int jk = 0; jk < nucl_num * (cord_num + 1); jk++) {
					for (int ml = 0; ml < 4 * elec_num; ml++) {

						// Single reduction
						int index_n_base = jk * stride_k_n + nw * stride_nw_n;
						int index_e_base =
							ml + i * stride_i_e + nw * stride_nw_e;
						sum = 0.;
#pragma omp simd reduction(+ : sum)
						for (int n = 0; n < elec_num; n++) {
							sum += een_rescaled_e_deriv_e[index_e_base +
														  n * stride_n_e] *
								   een_rescaled_n[index_n_base + n];
						}
						dtmp_c[ml + jk * stride_k_d + i * stride_i_d +
							   nw * stride_nw_d] = sum;
					}
				}
			}
		}
	}

	return QMCKL_SUCCESS;
}

//**********
// PROVIDE
//**********

qmckl_exit_code qmckl_provide_tmp_c_omp_offload(qmckl_context context) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	/* Check if dim_cord_vect is provided */
	qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.tmp_c_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			free(ctx->jastrow.tmp_c);
			ctx->jastrow.tmp_c = NULL;
		}

		/* Allocate array */
		if (ctx->jastrow.tmp_c == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size = (ctx->jastrow.cord_num) *
							(ctx->jastrow.cord_num + 1) * ctx->electron.num *
							ctx->nucleus.num * ctx->electron.walker.num *
							sizeof(double);
			double *tmp_c = (double *)qmckl_malloc(context, mem_info);

			if (tmp_c == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_provide_tmp_c_omp_offload", NULL);
			}
			ctx->jastrow.tmp_c = tmp_c;
		}

		rc = qmckl_compute_tmp_c_omp_offload(
			context, ctx->jastrow.cord_num, ctx->electron.num, ctx->nucleus.num,
			ctx->electron.walker.num, ctx->jastrow.een_rescaled_e,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.tmp_c);
	}

	ctx->jastrow.tmp_c_date = ctx->date;

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_dtmp_c_omp_offload(qmckl_context context) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	/* Check if dim_cord_vect is provided */
	qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.dtmp_c_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			free(ctx->jastrow.dtmp_c);
			ctx->jastrow.dtmp_c = NULL;
		}

		/* Allocate array */
		if (ctx->jastrow.dtmp_c == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size = (ctx->jastrow.cord_num) *
							(ctx->jastrow.cord_num + 1) * 4 *
							ctx->electron.num * ctx->nucleus.num *
							ctx->electron.walker.num * sizeof(double);
			double *dtmp_c = (double *)qmckl_malloc(context, mem_info);

			if (dtmp_c == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_provide_dtmp_c_omp_offload", NULL);
			}
			ctx->jastrow.dtmp_c = dtmp_c;
		}

		rc = qmckl_compute_dtmp_c_omp_offload(
			context, ctx->jastrow.cord_num, ctx->electron.num, ctx->nucleus.num,
			ctx->electron.walker.num, ctx->jastrow.een_rescaled_e_deriv_e,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.dtmp_c);

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->jastrow.dtmp_c_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_provide_factor_een_deriv_e_omp_offload(qmckl_context context) {

	qmckl_exit_code rc;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_een_rescaled_e(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_een_rescaled_n(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if en rescaled distance is provided */
	rc = qmckl_provide_een_rescaled_e_deriv_e(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_een_rescaled_n_deriv_e(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_cord_vect_full(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if en rescaled distance derivatives is provided */
	rc = qmckl_provide_lkpm_combined_index(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if tmp_c is provided */
	rc = qmckl_provide_tmp_c_omp_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Check if dtmp_c is provided */
	rc = qmckl_provide_dtmp_c_omp_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	/* Compute if necessary */
	if (ctx->date > ctx->jastrow.factor_een_deriv_e_date) {

		if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
			free(ctx->jastrow.factor_een_deriv_e);
			ctx->jastrow.factor_een_deriv_e = NULL;
		}

		/* Allocate array */
		if (ctx->jastrow.factor_een_deriv_e == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num *
							sizeof(double);
			double *factor_een_deriv_e =
				(double *)qmckl_malloc(context, mem_info);

			if (factor_een_deriv_e == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_provide_factor_"
									  "een_deriv_e_omp_offload",
									  NULL);
			}
			ctx->jastrow.factor_een_deriv_e = factor_een_deriv_e;
		}

		rc = qmckl_compute_factor_een_deriv_e(
			context, ctx->electron.walker.num, ctx->electron.num,
			ctx->nucleus.num, ctx->jastrow.cord_num, ctx->jastrow.dim_cord_vect,
			ctx->jastrow.cord_vect_full, ctx->jastrow.lkpm_combined_index,
			ctx->jastrow.tmp_c, ctx->jastrow.dtmp_c,
			ctx->jastrow.een_rescaled_n, ctx->jastrow.een_rescaled_n_deriv_e,
			ctx->jastrow.factor_een_deriv_e);
		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->jastrow.factor_een_deriv_e_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

//**********
// GET
//**********

qmckl_exit_code qmckl_get_jastrow_tmp_c_omp_offload(qmckl_context context,
													double *const tmp_c) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_dim_cord_vect(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	rc = qmckl_provide_cord_vect_full(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	rc = qmckl_provide_tmp_c_omp_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1) *
				 ctx->electron.num * ctx->nucleus.num *
				 ctx->electron.walker.num;
	memcpy(tmp_c, ctx->jastrow.tmp_c, sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_dtmp_c_omp_offload(qmckl_context context,
													 double *const dtmp_c) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_dim_cord_vect(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	rc = qmckl_provide_cord_vect_full(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	rc = qmckl_provide_dtmp_c_omp_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1) * 4 *
				 ctx->electron.num * ctx->nucleus.num *
				 ctx->electron.walker.num;
	memcpy(dtmp_c, ctx->jastrow.dtmp_c, sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_factor_een_deriv_e_omp_offload(
	qmckl_context context, double *const factor_een_deriv_e,
	const int64_t size_max) {
	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return QMCKL_NULL_CONTEXT;
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_factor_een_deriv_e_omp_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
	if (size_max < sze) {
		return qmckl_failwith(
			context, QMCKL_INVALID_ARG_3,
			"qmckl_get_jastrow_factor_een_deriv_e_omp_offload",
			"Array too small. Expected 4*walk_num*elec_num");
	}
	memcpy(factor_een_deriv_e, ctx->jastrow.factor_een_deriv_e,
		   sze * sizeof(double));

	return QMCKL_SUCCESS;
}
