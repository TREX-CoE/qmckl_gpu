#include "include/qmckl_ao.h"
#define qmckl_ten3(t, i, j, k) t.data[(i) + t.size[0]*((j) + t.size[1]*(k))]

//**********
// COMPUTE
//**********

#pragma acc routine seq
qmckl_exit_code qmckl_ao_polynomial_transp_vgl_hpc_acc_offload(
	const qmckl_context context, const double *restrict X,
	const double *restrict R, const int32_t lmax, int64_t *restrict n,
	const int64_t ldl, double *restrict const VGL, const int64_t ldv) {
	if (lmax < 0)
		return QMCKL_INVALID_ARG_4;
	if (ldl < 3)
		return QMCKL_INVALID_ARG_7;

	int32_t m;

	const int32_t size_max = (lmax + 1) * (lmax + 2) * (lmax + 3) / 6;
	if (ldv < size_max)
		return QMCKL_INVALID_ARG_9;

	double *restrict const vgl1 = VGL;
	double *restrict const vgl2 = VGL + ldv;
	double *restrict const vgl3 = VGL + (ldv << 1);
	double *restrict const vgl4 = VGL + ldv + (ldv << 1);
	double *restrict const vgl5 = VGL + (ldv << 2);

	const double Y[3] = {X[0] - R[0], X[1] - R[1], X[2] - R[2]};

	for (int32_t k = 0; k < 4; ++k) {
		vgl2[k] = 0.0;
		vgl3[k] = 0.0;
		vgl4[k] = 0.0;
		vgl5[k] = 0.0;
	}
	vgl1[0] = 1.0;
	vgl1[1] = Y[0];
	vgl1[2] = Y[1];
	vgl1[3] = Y[2];

	vgl2[1] = 1.0;
	vgl3[2] = 1.0;
	vgl4[3] = 1.0;
	m = 4;

	double pow_1, pow_2, pow_3;
	double dd = 2.0;

	for (int32_t d = 2; d <= lmax; ++d) {
		double da = dd;

		for (int32_t a = d; a >= 0; --a) {
			double db = dd - da;

			for (int32_t b = d - a; b >= 0; --b) {
				const int32_t c = d - a - b;
				const double dc = dd - da - db;

				// Compute pow_1 up to the power of a (a loop)
				pow_1 = 1.0f;
				for (int i = 0; i < a; ++i) {
					pow_1 = pow_1 * Y[0];
				}
				// Compute pow_2 up to the power of b (b loop)
				pow_2 = 1.0f;
				for (int i = 0; i < b; ++i) {
					pow_2 = pow_2 * Y[1];
				}
				// Compute pow_3 up to the power of c (c loop)
				pow_3 = 1.0f;
				for (int i = 0; i < c; ++i) {
					pow_3 = pow_3 * Y[2];
				}

				double xy = pow_1 * pow_2;
				double yz = pow_2 * pow_3;
				double xz = pow_1 * pow_3;

				vgl1[m] = xy * pow_3;

				xy *= dc;
				xz *= db;
				yz *= da;

				if (a >= 2)
					vgl2[m] = pow_1 / Y[0] * yz;
				else
					vgl2[m] = 1.0 * yz;
				if (b >= 2)
					vgl3[m] = pow_2 / Y[1] * xz;
				else
					vgl3[m] = 1.0 * xz;
				if (c >= 2)
					vgl4[m] = pow_3 / Y[2] * xy;
				else
					vgl4[m] = 1.0 * xy;

				double pow_1_tmp = (a >= 3) ? pow_1 : 1.;
				double pow_2_tmp = (b >= 3) ? pow_2 : 1.;
				double pow_3_tmp = (c >= 3) ? pow_3 : 1.;

				vgl5[m] = (da - 1.) * pow_1_tmp * yz +
						  (db - 1.) * pow_2_tmp * xz +
						  (dc - 1.) * pow_3_tmp * xy;

				db -= 1.0;
				++m;
			}
			da -= 1.0;
		}
		dd += 1.0;
	}

	*n = m;
	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_ao_vgl_gaussian_acc_offload(
	const qmckl_context context, const int64_t ao_num, const int64_t shell_num,
	const int32_t *restrict prim_num_per_nucleus, const int64_t point_num,
	const int64_t nucl_num, const double *restrict coord,
	const double *restrict nucl_coord, const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	const qmckl_matrix expo_per_nucleus, const qmckl_tensor coef_per_nucleus,
	double *restrict const ao_vgl) {
	int32_t lstart_h[32];
	for (int32_t l = 0; l < 32; ++l) {
		lstart_h[l] = l * (l + 1) * (l + 2) / 6;
	}

	int64_t *ao_index = malloc((shell_num + 1) * sizeof(int64_t));

	int64_t size_max = 0;
	int64_t prim_max = 0;
	int64_t shell_max = 0;
	int64_t k = 0;
	for (int inucl = 0; inucl < nucl_num; ++inucl) {
		prim_max = prim_num_per_nucleus[inucl] > prim_max
					   ? prim_num_per_nucleus[inucl]
					   : prim_max;
		shell_max = nucleus_shell_num[inucl] > shell_max
						? nucleus_shell_num[inucl]
						: shell_max;
		const int64_t ishell_start = nucleus_index[inucl];
		const int64_t ishell_end =
			nucleus_index[inucl] + nucleus_shell_num[inucl];
		for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {
			const int l = shell_ang_mom[ishell];
			ao_index[ishell] = k;
			k += lstart_h[l + 1] - lstart_h[l];
			size_max = size_max < lstart_h[l + 1] ? lstart_h[l + 1] : size_max;
		}
	}
	ao_index[shell_num] = ao_num + 1;

	/* Don't compute polynomials when the radial part is zero. */
	double cutoff = -log(1.e-12);

	double *poly_vgl_shared =
		(double *)malloc(point_num * 5 * size_max * sizeof(double));

	double *coef_mat =
		(double *)malloc(nucl_num * shell_max * prim_max * sizeof(double));
	for (int i = 0; i < nucl_num; ++i) {
		for (int j = 0; j < shell_max; ++j) {
			for (int k = 0; k < prim_max; ++k) {
				coef_mat[i * shell_max * prim_max + j * prim_max + k] =
					qmckl_ten3(coef_per_nucleus, k, j, i);
			}
		}
	}

	// WARNING This probably breaks the restrict on expo_per_nucleus.data
	double *expo_per_nucleus_data = expo_per_nucleus.data;
	int expo_per_nucleus_size_0 = expo_per_nucleus.size[0];
	int expo_per_nucleus_size_1 = expo_per_nucleus.size[1];

#pragma acc data copyin(                                                       \
		prim_num_per_nucleus[0 : nucl_num], coord[0 : 3 * point_num],          \
			nucl_coord[0 : 3 * nucl_num], nucleus_index[0 : nucl_num],         \
			nucleus_shell_num[0 : nucl_num], nucleus_range[0 : nucl_num],      \
			nucleus_max_ang_mom[0 : nucl_num], shell_ang_mom[0 : shell_num],   \
			ao_factor[0 : ao_num],                                             \
			expo_per_nucleus_data[0 : expo_per_nucleus_size_0 *                \
									  expo_per_nucleus_size_1],                \
			coef_mat[0 : nucl_num * shell_max * prim_max],                     \
			ao_index[0 : shell_num + 1])                                       \
	create(poly_vgl_shared[0 : point_num * 5 * size_max])
	{

#pragma acc parallel loop independent gang worker vector
		for (int64_t ipoint = 0; ipoint < point_num; ++ipoint) {

			double *poly_vgl = &(poly_vgl_shared[ipoint * 5 * size_max]);

			int32_t lstart[32];
			for (int32_t l = 0; l < 32; ++l) {
				lstart[l] = l * (l + 1) * (l + 2) / 6;
			}
			double poly_vgl_l1[4][4] = {{1.0, 0.0, 0.0, 0.0},
										{0.0, 1.0, 0.0, 0.0},
										{0.0, 0.0, 1.0, 0.0},
										{0.0, 0.0, 0.0, 1.0}};
			double poly_vgl_l2[5][10] = {
				{1., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
				{0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
				{0., 0., 1., 0., 0., 0., 0., 0., 0., 0.},
				{0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
				{0., 0., 0., 0., 2., 0., 0., 2., 0., 2.}};

			const double e_coord[3] = {coord[ipoint], coord[ipoint + point_num],
									   coord[ipoint + 2 * point_num]};

			for (int64_t inucl = 0; inucl < nucl_num; ++inucl) {
				const double n_coord[3] = {nucl_coord[inucl],
										   nucl_coord[inucl + nucl_num],
										   nucl_coord[inucl + 2 * nucl_num]};

				/* Test if the point is in the range of the
				 * nucleus */
				const double x = e_coord[0] - n_coord[0];
				const double y = e_coord[1] - n_coord[1];
				const double z = e_coord[2] - n_coord[2];

				const double r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				int64_t n_poly;
				switch (nucleus_max_ang_mom[inucl]) {
				case 0:
					break;

				case 1:
					poly_vgl_l1[0][1] = x;
					poly_vgl_l1[0][2] = y;
					poly_vgl_l1[0][3] = z;
					break;

				case 2:
					poly_vgl_l2[0][1] = x;
					poly_vgl_l2[0][2] = y;
					poly_vgl_l2[0][3] = z;
					poly_vgl_l2[0][4] = x * x;
					poly_vgl_l2[0][5] = x * y;
					poly_vgl_l2[0][6] = x * z;
					poly_vgl_l2[0][7] = y * y;
					poly_vgl_l2[0][8] = y * z;
					poly_vgl_l2[0][9] = z * z;
					poly_vgl_l2[1][4] = x + x;
					poly_vgl_l2[1][5] = y;
					poly_vgl_l2[1][6] = z;
					poly_vgl_l2[2][5] = x;
					poly_vgl_l2[2][7] = y + y;
					poly_vgl_l2[2][8] = z;
					poly_vgl_l2[3][6] = x;
					poly_vgl_l2[3][8] = y;
					poly_vgl_l2[3][9] = z + z;
					break;

				default:
					qmckl_ao_polynomial_transp_vgl_hpc_acc_offload(
						context, e_coord, n_coord, nucleus_max_ang_mom[inucl],
						&n_poly, (int64_t)3, poly_vgl, size_max);

					break;
				}

				/* Compute all exponents */

				int64_t nidx = 0;
				int base_idx = inucl * expo_per_nucleus_size_0;
				for (int64_t iprim = 0; iprim < prim_num_per_nucleus[inucl];
					 ++iprim) {
					const double v =
						expo_per_nucleus_data[base_idx + iprim] * r2;
					if (v <= cutoff) {
						++nidx;
					} else {
						break;
					}
				}

				const int64_t ishell_start = nucleus_index[inucl];
				const int64_t ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl];

				double ce_mat_0 = 0.;
				double ce_mat_1 = 0.;
				double ce_mat_2 = 0.;
				double ce_mat_3 = 0.;
				double ce_mat_4 = 0.;

				double exp_mat_0 = 0.;
				double exp_mat_1 = 0.;
				double exp_mat_2 = 0.;
				double exp_mat_3 = 0.;
				double exp_mat_4 = 0.;

				for (int64_t ishell = ishell_start; ishell < ishell_end;
					 ++ishell) {
					ce_mat_0 = 0.;
					ce_mat_1 = 0.;
					ce_mat_2 = 0.;
					ce_mat_3 = 0.;
					ce_mat_4 = 0.;

					for (int k = 0; k < nidx; ++k) {
						exp_mat_0 =
							exp(-(expo_per_nucleus_data[base_idx + k] * r2));
						double f =
							expo_per_nucleus_data[base_idx + k] * exp_mat_0;
						f = -f - f;
						exp_mat_1 = f * x;
						exp_mat_2 = f * y;
						exp_mat_3 = f * z;
						exp_mat_4 =
							f *
							(3.0 -
							 2.0 * (expo_per_nucleus_data[base_idx + k] * r2));
						if (coef_mat[inucl * shell_max * prim_max +
									 (ishell - ishell_start) * prim_max + k] !=
							0.) {
							ce_mat_0 +=
								coef_mat[inucl * shell_max * prim_max +
										 (ishell - ishell_start) * prim_max +
										 k] *
								exp_mat_0;
							ce_mat_1 +=
								coef_mat[inucl * shell_max * prim_max +
										 (ishell - ishell_start) * prim_max +
										 k] *
								exp_mat_1;
							ce_mat_2 +=
								coef_mat[inucl * shell_max * prim_max +
										 (ishell - ishell_start) * prim_max +
										 k] *
								exp_mat_2;
							ce_mat_3 +=
								coef_mat[inucl * shell_max * prim_max +
										 (ishell - ishell_start) * prim_max +
										 k] *
								exp_mat_3;
							ce_mat_4 +=
								coef_mat[inucl * shell_max * prim_max +
										 (ishell - ishell_start) * prim_max +
										 k] *
								exp_mat_4;
						}
					}

					const double s1 = ce_mat_0;
					if (s1 == 0.0)
						continue;
					const double s2 = ce_mat_1;
					const double s3 = ce_mat_2;
					const double s4 = ce_mat_3;
					const double s5 = ce_mat_4;

					const int64_t k = ao_index[ishell];
					double *restrict const ao_vgl_1 =
						ao_vgl + ipoint * 5 * ao_num + k;

					const int32_t l = shell_ang_mom[ishell];
					const int32_t n = lstart[l + 1] - lstart[l];

					double *restrict const ao_vgl_2 = ao_vgl_1 + ao_num;
					double *restrict const ao_vgl_3 = ao_vgl_1 + (ao_num << 1);
					double *restrict const ao_vgl_4 =
						ao_vgl_1 + (ao_num << 1) + ao_num;
					double *restrict const ao_vgl_5 = ao_vgl_1 + (ao_num << 2);

					double *restrict poly_vgl_1 = NULL;
					double *restrict poly_vgl_2 = NULL;
					double *restrict poly_vgl_3 = NULL;
					double *restrict poly_vgl_4 = NULL;
					double *restrict poly_vgl_5 = NULL;
					if (nidx > 0) {
						const double *restrict f = ao_factor + k;
						const int64_t idx = lstart[l];

						switch (nucleus_max_ang_mom[inucl]) {
						case 0:
							break;
						case 1:
							poly_vgl_1 = &(poly_vgl_l1[0][idx]);
							poly_vgl_2 = &(poly_vgl_l1[1][idx]);
							poly_vgl_3 = &(poly_vgl_l1[2][idx]);
							poly_vgl_4 = &(poly_vgl_l1[3][idx]);
							break;
						case 2:
							poly_vgl_1 = &(poly_vgl_l2[0][idx]);
							poly_vgl_2 = &(poly_vgl_l2[1][idx]);
							poly_vgl_3 = &(poly_vgl_l2[2][idx]);
							poly_vgl_4 = &(poly_vgl_l2[3][idx]);
							poly_vgl_5 = &(poly_vgl_l2[4][idx]);
							break;
						default:
							poly_vgl_1 = &(poly_vgl[idx]);
							poly_vgl_2 = &(poly_vgl[size_max + idx]);
							poly_vgl_3 = &(poly_vgl[2 * size_max + idx]);
							poly_vgl_4 = &(poly_vgl[3 * size_max + idx]);
							poly_vgl_5 = &(poly_vgl[4 * size_max + idx]);
						}
						switch (n) {
						case (1):
							ao_vgl_1[0] = s1 * f[0];
							ao_vgl_2[0] = s2 * f[0];
							ao_vgl_3[0] = s3 * f[0];
							ao_vgl_4[0] = s4 * f[0];
							ao_vgl_5[0] = s5;
							break;
						case (3):
							for (int il = 0; il < 3; ++il) {
								ao_vgl_1[il] = poly_vgl_1[il] * s1 * f[il];
								ao_vgl_2[il] = (poly_vgl_2[il] * s1 +
												poly_vgl_1[il] * s2) *
											   f[il];
								ao_vgl_3[il] = (poly_vgl_3[il] * s1 +
												poly_vgl_1[il] * s3) *
											   f[il];
								ao_vgl_4[il] = (poly_vgl_4[il] * s1 +
												poly_vgl_1[il] * s4) *
											   f[il];
								ao_vgl_5[il] = (poly_vgl_1[il] * s5 +
												2.0 * (poly_vgl_2[il] * s2 +
													   poly_vgl_3[il] * s3 +
													   poly_vgl_4[il] * s4)) *
											   f[il];
							}
							break;
						case (6):
							for (int il = 0; il < 6; ++il) {
								ao_vgl_1[il] = poly_vgl_1[il] * s1 * f[il];
								ao_vgl_2[il] = (poly_vgl_2[il] * s1 +
												poly_vgl_1[il] * s2) *
											   f[il];
								ao_vgl_3[il] = (poly_vgl_3[il] * s1 +
												poly_vgl_1[il] * s3) *
											   f[il];
								ao_vgl_4[il] = (poly_vgl_4[il] * s1 +
												poly_vgl_1[il] * s4) *
											   f[il];
								ao_vgl_5[il] =
									(poly_vgl_5[il] * s1 + poly_vgl_1[il] * s5 +
									 2.0 * (poly_vgl_2[il] * s2 +
											poly_vgl_3[il] * s3 +
											poly_vgl_4[il] * s4)) *
									f[il];
							}
							break;
						default:
							for (int il = 0; il < n; ++il) {
								ao_vgl_1[il] = poly_vgl_1[il] * s1 * f[il];
								ao_vgl_2[il] = (poly_vgl_2[il] * s1 +
												poly_vgl_1[il] * s2) *
											   f[il];
								ao_vgl_3[il] = (poly_vgl_3[il] * s1 +
												poly_vgl_1[il] * s3) *
											   f[il];
								ao_vgl_4[il] = (poly_vgl_4[il] * s1 +
												poly_vgl_1[il] * s4) *
											   f[il];
								ao_vgl_5[il] =
									(poly_vgl_5[il] * s1 + poly_vgl_1[il] * s5 +
									 2.0 * (poly_vgl_2[il] * s2 +
											poly_vgl_3[il] * s3 +
											poly_vgl_4[il] * s4)) *
									f[il];
							}
							break;
						}
					} else {
						for (int64_t il = 0; il < n; ++il) {
							ao_vgl_1[il] = 0.0;
							ao_vgl_2[il] = 0.0;
							ao_vgl_3[il] = 0.0;
							ao_vgl_4[il] = 0.0;
							ao_vgl_5[il] = 0.0;
						}
					}
				}
			}
		}
	}

	free(ao_index);
	free(coef_mat);
	free(poly_vgl_shared);

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ao_vgl_acc_offload(qmckl_context context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_vgl", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED, "qmckl_ao_vgl",
							  NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

		qmckl_exit_code rc;

		/* Allocate array */
		if (ctx->ao_basis.ao_vgl == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double);
			double *ao_vgl = (double *)qmckl_malloc(context, mem_info);

			if (ao_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_ao_basis_ao_vgl", NULL);
			}
			ctx->ao_basis.ao_vgl = ao_vgl;
		}

#pragma acc enter data copyin(ctx)
#pragma acc enter data create(                                                 \
		ctx->ao_basis.ao_vgl[0 : ctx->point.num * 5 * ctx->ao_basis.ao_num])

		if (ctx->ao_basis.type == 'G') {
			rc = qmckl_compute_ao_vgl_gaussian_acc_offload(
				context, ctx->ao_basis.ao_num, ctx->ao_basis.shell_num,
				ctx->ao_basis.prim_num_per_nucleus, ctx->point.num,
				ctx->nucleus.num, ctx->point.coord.data,
				ctx->nucleus.coord.data, ctx->ao_basis.nucleus_index,
				ctx->ao_basis.nucleus_shell_num, ctx->ao_basis.nucleus_range,
				ctx->ao_basis.nucleus_max_ang_mom, ctx->ao_basis.shell_ang_mom,
				ctx->ao_basis.ao_factor, ctx->ao_basis.expo_per_nucleus,
				ctx->ao_basis.coef_per_nucleus, ctx->ao_basis.ao_vgl);
		} else {
			printf("ERROR: basis type different from 'G' are not "
				   "supported in this "
				   "function's version\n");
			rc = QMCKL_FAILURE;
		}

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->ao_basis.ao_vgl_date = ctx->date;
		// TG: uncomment ao_vgl are needed on the host before the computation of
		// mo_vgl
		// #pragma acc update host(ctx->ao_basis.ao_vgl [0:ctx->point.num * 5 *
		// ctx->ao_basis.ao_num])
	}

	return QMCKL_SUCCESS;
}

//**********
// GET
//**********

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_acc_offload(qmckl_context context,
													  double *const ao_vgl,
													  const int64_t size_max) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_vgl", NULL);
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_ao_vgl_acc_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_ao_vgl",
							  "input array too small");
	}
	memcpy(ao_vgl, ctx->ao_basis.ao_vgl, (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_inplace_acc_offload(
	qmckl_context context, double *const ao_vgl, const int64_t size_max) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_vgl", NULL);
	}

	qmckl_exit_code rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_ao_vgl",
							  "input array too small");
	}

	rc = qmckl_context_touch(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	double *old_array = ctx->ao_basis.ao_vgl;

	ctx->ao_basis.ao_vgl = ao_vgl;

	rc = qmckl_provide_ao_vgl_acc_offload(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	ctx->ao_basis.ao_vgl = old_array;

	return QMCKL_SUCCESS;
}
