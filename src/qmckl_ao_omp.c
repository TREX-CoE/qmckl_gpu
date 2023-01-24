#include "include/qmckl_ao.h"

//**********
// COMPUTE
//**********

/* shell_vgl */

qmckl_exit_code qmckl_compute_ao_basis_shell_gaussian_vgl_device(
	qmckl_context_device context, int prim_num, int shell_num, int point_num,
	int nucl_num, int64_t *nucleus_shell_num, int64_t *nucleus_index,
	double *nucleus_range, int64_t *shell_prim_index, int64_t *shell_prim_num,
	double *coord, double *nucl_coord, double *expo, double *coef_normalized,
	double *shell_vgl) {

	int ishell_start, ishell_end;
	int iprim_start, iprim_end;
	double x, y, z, two_a, ar2, r2, v, cutoff;

	qmckl_exit_code info = QMCKL_SUCCESS;

	// Don't compute exponentials when the result will be almost zero.
	// TODO : Use numerical precision here
	cutoff = -log(1e-12);

#pragma omp target is_device_ptr(                                              \
	nucleus_shell_num, nucleus_index, nucleus_range, shell_prim_index,         \
	shell_prim_num, coord, nucl_coord, expo, coef_normalized, shell_vgl)
	{

#pragma omp teams distribute parallel for simd collapse(2)
		for (int ipoint = 0; ipoint < point_num; ipoint++) {

			for (int inucl = 0; inucl < nucl_num; inucl++) {

				x = coord[ipoint] - nucl_coord[inucl];
				y = coord[ipoint + point_num] - nucl_coord[inucl + nucl_num];
				z = coord[ipoint + 2 * point_num] -
					nucl_coord[inucl + 2 * nucl_num];

				r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				// C is zero-based, so shift bounds by one
				ishell_start = nucleus_index[inucl] + 1;
				ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl];

				for (int ishell = ishell_start; ishell < ishell_end; ishell++) {

					shell_vgl[ishell, 1, ipoint] = 0;
					shell_vgl[ishell, 2, ipoint] = 0;
					shell_vgl[ishell, 3, ipoint] = 0;
					shell_vgl[ishell, 4, ipoint] = 0;
					shell_vgl[ishell, 5, ipoint] = 0;

					iprim_start = shell_prim_index[ishell] + 1;
					iprim_end =
						shell_prim_index[ishell] + shell_prim_num[ishell];

					for (int iprim = iprim_start; iprim < iprim_end; iprim++) {

						ar2 = expo[iprim] * r2;
						if (ar2 > cutoff) {
							continue;
						}

						v = coef_normalized[iprim] * exp(-ar2);
						two_a = -2 * expo[iprim] * v;

						shell_vgl[ishell + 0 * shell_num +
								  ipoint * shell_num * 5] =
							shell_vgl[ishell + 0 * shell_num +
									  ipoint * shell_num * 5] +
							v;

						shell_vgl[ishell + 1 * shell_num +
								  ipoint * shell_num * 5] =
							shell_vgl[ishell + 1 * shell_num +
									  ipoint * shell_num * 5] +
							two_a * x;

						shell_vgl[ishell + 2 * shell_num +
								  ipoint * shell_num * 5] =
							shell_vgl[ishell + 2 * shell_num +
									  ipoint * shell_num * 5] +
							two_a * y;

						shell_vgl[ishell + 3 * shell_num +
								  ipoint * shell_num * 5] =
							shell_vgl[ishell + 3 * shell_num +
									  ipoint * shell_num * 5] +
							two_a * z;

						shell_vgl[ishell + 4 * shell_num +
								  ipoint * shell_num * 5] =
							shell_vgl[ishell + 4 * shell_num +
									  ipoint * shell_num * 5] +
							two_a * (3 - 2 * ar2);
					}
				}
			}
		}
	}

	return info;
}

/* ao_vgl_gaussian */

qmckl_exit_code qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context context, const int64_t ao_num, const int64_t shell_num,
	const int32_t *restrict prim_num_per_nucleus, const int64_t point_num,
	const int64_t nucl_num, const double *restrict coord,
	const double *restrict nucl_coord, const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_vgl) {

	double *e_coord, *n_coord;
	int64_t n_poly;
	int64_t l, il, k;
	int64_t ipoint, inucl, ishell;
	int64_t ishell_start, ishell_end;
	int64_t *lstart;
	double x, y, z, r2;
	double cutoff;
	int64_t qmckl_ao_polynomial_vgl_doc_f;
	int64_t size_max = 0;

	double *poly_vgl;
	int64_t *powers;
	int64_t *ao_index;

	qmckl_exit_code rc;
	int lmax, c;

	int device_id = qmckl_get_device_id(context);
	qmckl_memory_info_struct info;

	info.size = sizeof(int64_t) * 21;
	lstart = qmckl_malloc_device(context, info);

	info.size = sizeof(double) * 8;
	e_coord = qmckl_malloc_device(context, info);
	info.size = sizeof(double) * 8;
	n_coord = qmckl_malloc_device(context, info);

	info.size = sizeof(double) * 5 * ao_num;
	poly_vgl = qmckl_malloc_device(context, info);
	info.size = sizeof(int64_t) * 3 * ao_num;
	powers = qmckl_malloc_device(context, info);
	info.size = sizeof(int64_t) * ao_num;
	ao_index = qmckl_malloc_device(context, info);

	// Specific calling function
	lmax = -1;
#pragma is_device_ptr(nucleus_max_ang_mom)
	{
#pragma omp target
		{
			for (int i = 0; i < nucl_num; i++) {
				if (lmax < nucleus_max_ang_mom[i]) {
					lmax = nucleus_max_ang_mom[i];
				}
			}
		}
	}
	info.size = 3 * (lmax + 2) * sizeof(double);
	double *pows = qmckl_malloc_device(context, info);

#pragma omp target is_device_ptr(lstart)
	{
		for (l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	k = 0;
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
								 shell_ang_mom, ao_index, lstart)
	{
		for (inucl = 0; inucl < nucl_num; inucl++) {
			ishell_start = nucleus_index[inucl] + 1;
			ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl];
			for (ishell = ishell_start; ishell < ishell_end; ishell++) {
				l = shell_ang_mom[ishell];
				ao_index[ishell] = k;
				k = k + lstart[l + 1] - lstart[l];
			}
		}
	}

#pragma omp target is_device_ptr(                                              \
	ao_vgl, lstart, e_coord, n_coord, ao_index, ao_factor, coord,              \
	nucleus_max_ang_mom, nucleus_index, nucleus_shell_num, shell_vgl,          \
	poly_vgl, nucl_coord, pows, powers, shell_ang_mom)
	{
#pragma omp teams distribute parallel for
		for (ipoint = 0; ipoint < point_num; ipoint++) {

			e_coord[0] = coord[0 * point_num + ipoint];
			e_coord[1] = coord[1 * point_num + ipoint];
			e_coord[2] = coord[2 * point_num + ipoint];

			for (inucl = 0; inucl < nucl_num; inucl++) {

				n_coord[0] = nucl_coord[0 * nucl_num + inucl];
				n_coord[1] = nucl_coord[1 * nucl_num + inucl];
				n_coord[2] = nucl_coord[2 * nucl_num + inucl];

				x = e_coord[0] - n_coord[0];
				y = e_coord[1] - n_coord[1];
				z = e_coord[2] - n_coord[2];

				r2 = x * x + y * y + z * z;

				// Beginning of ao_polynomial computation (now inlined)
				double Y1, Y2, Y3;
				double xy, yz, xz;
				int c;
				double da, db, dc, dd;

				Y1 = e_coord[0] - n_coord[0];
				Y2 = e_coord[1] - n_coord[1];
				Y3 = e_coord[2] - n_coord[2];

				if (nucleus_max_ang_mom[inucl] == 0) {

					poly_vgl[0] = 1;
					poly_vgl[1 * 3 + 0] = 0;
					poly_vgl[2 * 3 + 0] = 0;
					poly_vgl[3 * 3 + 0] = 0;
					poly_vgl[4 * 3 + 0] = 0;

					for (int i = 0; i < 3; i++) {
						powers[0 * 3 + i] = 0;
					}

					n_poly = 0;

				}

				else if (nucleus_max_ang_mom[inucl] > 0) {
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							pows[i * (nucleus_max_ang_mom[inucl] + 2) + j] = 1;
						}
					}

					/**/
					for (int i = 3; i < nucleus_max_ang_mom[inucl] + 2; i++) {
						pows[0 * (nucleus_max_ang_mom[inucl] + 2) + i] =
							pows[0 * (nucleus_max_ang_mom[inucl] + 2) +
								 (i - 1)] *
							Y1;
						pows[1 * (nucleus_max_ang_mom[inucl] + 2) + i] =
							pows[1 * (nucleus_max_ang_mom[inucl] + 2) +
								 (i - 1)] *
							Y2;
						pows[2 * (nucleus_max_ang_mom[inucl] + 2) + i] =
							pows[2 * (nucleus_max_ang_mom[inucl] + 2) +
								 (i - 1)] *
							Y3;
					}
					/**/

					for (int i = 0; i < 5; i++) {
						for (int j = 0; j < 4; j++) {
							poly_vgl[i * 3 + j] = 0;
						}
					}
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 4; j++) {
							powers[i * 3 + j] = 0;
						}
					}

					poly_vgl[0] = 1;

					powers[1] = 1;
					poly_vgl[1] = pows[0];
					poly_vgl[1 * 3 + 1] = 1;
					poly_vgl[2 * 3 + 2] = 1;

					n_poly = 4;
				}
				// l >= 2

				dd = 2;
				for (int d = 4; d < nucleus_max_ang_mom[inucl]; d++) {
					da = dd;
					for (int a = d; a > 1; a--) {
						db = dd - da;
						for (int b = d - a; b > 1; b--) {

							c = d - a - b;
							dc = dd - da - db;
							n_poly++;

							powers[n_poly * 3 + 0] = a;
							powers[n_poly * 3 + 1] = b;
							powers[n_poly * 3 + 2] = c;

							xy = pows[a] *
								 pows[1 * (nucleus_max_ang_mom[inucl] + 2) + b];
							yz =
								pows[1 * (nucleus_max_ang_mom[inucl] + 2) + b] *
								pows[2 * (nucleus_max_ang_mom[inucl] + 2) + c];
							yz = pows[a] *
								 pows[2 * (nucleus_max_ang_mom[inucl] + 2) + c];

							poly_vgl[n_poly * ao_num + 0] =
								xy *
								powers[3 * (nucleus_max_ang_mom[inucl] + 2) +
									   2];

							xy = dc * xy;
							xz = db * xz;
							yz = da * yz;

							poly_vgl[n_poly * 5 + 1] =
								pows[0 * (nucleus_max_ang_mom[inucl] + 2) + a -
									 1] *
								yz;
							poly_vgl[n_poly * 5 + 2] =
								pows[1 * (nucleus_max_ang_mom[inucl] + 2) +
									 (b - 1)] *
								xz;
							poly_vgl[n_poly * 5 + 3] =
								pows[2 * (nucleus_max_ang_mom[inucl] + 2) +
									 (c - 1)] *
								xy;

							poly_vgl[n_poly * ao_num + 4] =
								(da - 1) * pows[a - 2] * yz +
								(db - 1) *
									pows[1 * (nucleus_max_ang_mom[inucl] + 2) +
										 (b - 2)] *
									xz +
								(dc - 1) *
									pows[2 * (nucleus_max_ang_mom[inucl] + 2) +
										 (c - 2)] *
									xy;

							db = db - 1;
						}
						da = da - 1;
					}
					dd = dd + 1;
				}
				// End of ao_polynomial computation (now inlined)

				// Start of debug comment 2
				ishell_start = nucleus_index[inucl] + 1;
				ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl];

				for (ishell = ishell_start; ishell < ishell_end; ishell++) {
					k = ao_index[ishell];
					l = shell_ang_mom[ishell];
					for (il = lstart[l]; il < lstart[l + 1]; il++) {
						// value
						ao_vgl[ipoint * 5 * ao_num + 0 * ao_num + k] =
							poly_vgl[il * 5 + 0] *
							shell_vgl[ipoint * 5 * shell_num + 0 * shell_num +
									  ishell] *
							ao_factor[k];

						// Grad x
						ao_vgl[ipoint * 5 * ao_num + 1 * ao_num + k] =
							poly_vgl[il * 5 + 0] *
								shell_vgl[ipoint * 5 * shell_num +
										  0 * shell_num + ishell] +
							poly_vgl[il * 5 + 1] *
								shell_vgl[ipoint * 5 * shell_num +
										  1 * shell_num + ishell] *
								ao_factor[k];

						// grad y
						ao_vgl[ipoint * 5 * ao_num + 2 * ao_num + k] =
							poly_vgl[il * 5 + 2] *
								shell_vgl[ipoint * 5 * shell_num +
										  0 * shell_num + ishell] +
							poly_vgl[il * 5 + 1] *
								shell_vgl[ipoint * 5 * shell_num +
										  2 * shell_num + ishell] *
								ao_factor[k];

						// grad z
						ao_vgl[ipoint * 5 * ao_num + 3 * ao_num + k] =
							poly_vgl[il * 5 + 3] *
								shell_vgl[ipoint * 5 * shell_num +
										  0 * shell_num + ishell] +
							poly_vgl[il * 5 + 1] *
								shell_vgl[ipoint * 5 * shell_num +
										  3 * shell_num + ishell] *
								ao_factor[k];

						// Lapl_z

						ao_vgl[ipoint * 5 * ao_num + 4 * ao_num + k] =
							(poly_vgl[il * 5 + 4] *
								 shell_vgl[ipoint * 5 * shell_num +
										   0 * shell_num + ishell] +
							 poly_vgl[il * 5 + 0] *
								 shell_vgl[ipoint * 5 * shell_num +
										   0 * shell_num + ishell] +
							 2.0 * (poly_vgl[il * 5 + 1] *
										shell_vgl[ipoint * 5 * shell_num +
												  1 * shell_num + ishell] +
									poly_vgl[il * 5 + 2] *
										shell_vgl[ipoint * 5 * shell_num +
												  2 * shell_num + ishell] +
									poly_vgl[il * 5 + 3] *
										shell_vgl[ipoint * 5 * shell_num +
												  3 * shell_num + ishell])) *
							ao_factor[k];

						k = k + 1;
					}
				}
			}
		}
	}
	qmckl_free_device(context, lstart);
	qmckl_free_device(context, e_coord);
	qmckl_free_device(context, n_coord);
	qmckl_free_device(context, poly_vgl);
	qmckl_free_device(context, powers);
	qmckl_free_device(context, ao_index);

	qmckl_free_device(context, pows);

	return QMCKL_SUCCESS;
}

//**********
// PROVIDE
//**********

/* shell_vgl */

qmckl_exit_code
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
							  "qmckl_provide_ao_basis_shell_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.shell_vgl_date) {

		/* Allocate array */
		if (ctx->ao_basis.shell_vgl == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->ao_basis.shell_num * 5 * ctx->point.num * sizeof(double);
			double *shell_vgl =
				(double *)qmckl_malloc_device(context, mem_info);

			if (shell_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_ao_basis_shell_vgl_device", NULL);
			}
			ctx->ao_basis.shell_vgl = shell_vgl;
		}

		qmckl_exit_code rc;
		if (ctx->ao_basis.type == 'G') {
			rc = qmckl_compute_ao_basis_shell_gaussian_vgl_device(
				context, ctx->ao_basis.prim_num, ctx->ao_basis.shell_num,
				ctx->point.num, ctx->nucleus.num,
				ctx->ao_basis.nucleus_shell_num, ctx->ao_basis.nucleus_index,
				ctx->ao_basis.nucleus_range, ctx->ao_basis.shell_prim_index,
				ctx->ao_basis.shell_prim_num, ctx->point.coord.data,
				ctx->nucleus.coord.data, ctx->ao_basis.exponent,
				ctx->ao_basis.coefficient_normalized, ctx->ao_basis.shell_vgl);
		} else {
			return qmckl_failwith(context, QMCKL_FAILURE,
								  "compute_ao_basis_shell_vgl",
								  "Not yet implemented for basis type != 'G'");
		}
		if (rc != QMCKL_SUCCESS) {

			return rc;
		}

		ctx->ao_basis.shell_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

/* ao_vgl_gaussian */

qmckl_exit_code
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith((qmckl_context)context, QMCKL_NOT_PROVIDED,
							  "qmckl_ao_basis_ao_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

		qmckl_exit_code rc;

		/* Allocate array */
		if (ctx->ao_basis.ao_vgl == NULL) {

			qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
			mem_info.size =
				ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double);
			double *ao_vgl = (double *)qmckl_malloc_device(context, mem_info);

			if (ao_vgl == NULL) {
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
									  "qmckl_ao_basis_ao_vgl", NULL);
			}
			ctx->ao_basis.ao_vgl = ao_vgl;
		}

		/* Checking for shell_vgl */
		if (ctx->ao_basis.shell_vgl == NULL ||
			ctx->point.date > ctx->ao_basis.shell_vgl_date) {
			qmckl_provide_ao_basis_shell_vgl_device(context);
		}

		/* Compute ao_vgl_gaussian itself */
		if (ctx->ao_basis.type == 'G') {
			rc = qmckl_compute_ao_vgl_gaussian_device(
				context, ctx->ao_basis.ao_num, ctx->ao_basis.shell_num,
				ctx->ao_basis.prim_num_per_nucleus, ctx->point.num,
				ctx->nucleus.num, ctx->point.coord.data,
				ctx->nucleus.coord.data, ctx->ao_basis.nucleus_index,
				ctx->ao_basis.nucleus_shell_num, ctx->ao_basis.nucleus_range,
				ctx->ao_basis.nucleus_max_ang_mom, ctx->ao_basis.shell_ang_mom,
				ctx->ao_basis.ao_factor, ctx->ao_basis.shell_vgl,
				ctx->ao_basis.ao_vgl);
		} else {
			printf("Device pointers version of ao_vgl only "
				   "supports 'G' as its "
				   "ao_basis.type for now\n ");
		}

		if (rc != QMCKL_SUCCESS) {
			return rc;
		}

		ctx->ao_basis.ao_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS;
}

//**********
// GET
//**********

/* shell_vgl */

qmckl_exit_code
qmckl_get_ao_basis_shell_vgl_device(qmckl_context_device context,
									double *const shell_vgl,
									const int64_t size_max) {
	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_shell_vgl_device", NULL);
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_ao_basis_shell_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);
	// int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.shell_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_shell_vgl",
							  "input array too small");
	}
	qmckl_memcpy_D2D(context, shell_vgl, ctx->ao_basis.shell_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS;
}

/* ao_vgl_gaussian */

qmckl_exit_code qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
												 double *const ao_vgl,
												 const int64_t size_max) {

	if (qmckl_context_check((qmckl_context)context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_CONTEXT,
							  "qmckl_get_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_exit_code rc;

	rc = qmckl_provide_ao_basis_ao_vgl_device(context);
	if (rc != QMCKL_SUCCESS)
		return rc;

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith((qmckl_context)context, QMCKL_INVALID_ARG_3,
							  "qmckl_get_ao_basis_ao_vgl_device",
							  "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_vgl, ctx->ao_basis.ao_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS;
}
