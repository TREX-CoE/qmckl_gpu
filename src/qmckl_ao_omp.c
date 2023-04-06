#include "include/qmckl_ao.h"

//**********
// COMPUTE
//**********

/* shell_vgl */

qmckl_exit_code_device qmckl_compute_ao_basis_shell_gaussian_vgl_device(
	qmckl_context_device context, int prim_num, int shell_num, int point_num,
	int nucl_num, int64_t *nucleus_shell_num, int64_t *nucleus_index,
	double *nucleus_range, int64_t *shell_prim_index, int64_t *shell_prim_num,
	double *coord, double *nucl_coord, double *expo, double *coef_normalized,
	double *shell_vgl) {

	int ishell_start, ishell_end;
	int iprim_start, iprim_end;
	double x, y, z, two_a, ar2, r2, v, cutoff;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	// Don't compute exponentials when the result will be almost zero.
	// TODO : Use numerical precision here
	cutoff = 27.631021115928547; //-dlog(1.d-12)

#pragma omp target is_device_ptr(nucleus_shell_num, nucleus_index,             \
									 nucleus_range, shell_prim_index,          \
									 shell_prim_num, coord, nucl_coord, expo,  \
									 coef_normalized, shell_vgl)
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

				ishell_start = nucleus_index[inucl];
				ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;

				for (int ishell = ishell_start; ishell <= ishell_end;
					 ishell++) {

					shell_vgl[ishell + 0 * shell_num + ipoint * shell_num * 5] =
						0;
					shell_vgl[ishell + 1 * shell_num + ipoint * shell_num * 5] =
						0;
					shell_vgl[ishell + 2 * shell_num + ipoint * shell_num * 5] =
						0;
					shell_vgl[ishell + 3 * shell_num + ipoint * shell_num * 5] =
						0;
					shell_vgl[ishell + 4 * shell_num + ipoint * shell_num * 5] =
						0;

					iprim_start = shell_prim_index[ishell];
					iprim_end =
						shell_prim_index[ishell] + shell_prim_num[ishell] - 1;

					for (int iprim = iprim_start; iprim <= iprim_end; iprim++) {

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

/* ao_vgl */

qmckl_exit_code_device qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_vgl) {

	int64_t n_poly;
	int64_t *lstart;
	double cutoff = 27.631021115928547;

	double *poly_vgl_shared;
	int64_t *powers;
	int64_t *ao_index;

	qmckl_exit_code_device rc;

	lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Multiply "normal" size by point_num to affect subarrays to each thread
	poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * 5 * ao_num * point_num);

	ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);

	// Specific calling function
	int lmax = -1;
#pragma omp target is_device_ptr(nucleus_max_ang_mom) map(tofrom : lmax)
	{
		for (int i = 0; i < nucl_num; i++) {
			if (lmax < nucleus_max_ang_mom[i]) {
				lmax = nucleus_max_ang_mom[i];
			}
		}
	}
	// Multiply "normal" size by point_num to affect subarrays to each thread
	double *pows_shared = qmckl_malloc_device(
		context, sizeof(double) * (lmax + 3) * 3 * point_num);

#pragma omp target is_device_ptr(lstart)
	{
		for (int l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	int k = 1;
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
									 shell_ang_mom, ao_index, lstart)          \
	map(tofrom : k)
	{
		for (int inucl = 0; inucl < nucl_num; inucl++) {
			int ishell_start = nucleus_index[inucl];
			int ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
			for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
				int l = shell_ang_mom[ishell];
				ao_index[ishell] = k;
				k = k + lstart[l + 1] - lstart[l];
			}
		}
	}

#pragma omp target is_device_ptr(                                              \
		ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,       \
			nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared,      \
			nucl_coord, pows_shared, shell_ang_mom, nucleus_range)
	{
#pragma omp teams distribute parallel for simd // private(pows, poly_vgl)

		for (int ipoint = 0; ipoint < point_num; ipoint++) {

			// Compute addresses of subarrays from ipoint
			// This way, each thread can write to its own poly_vgl and pows
			// without any race condition
			double *poly_vgl = poly_vgl_shared + ipoint * 5 * ao_num;
			double *pows = pows_shared + ipoint * (lmax + 3) * 3;

			double e_coord_0 = coord[0 * point_num + ipoint];
			double e_coord_1 = coord[1 * point_num + ipoint];
			double e_coord_2 = coord[2 * point_num + ipoint];

			for (int inucl = 0; inucl < nucl_num; inucl++) {

				double n_coord_0 = nucl_coord[0 * nucl_num + inucl];
				double n_coord_1 = nucl_coord[1 * nucl_num + inucl];
				double n_coord_2 = nucl_coord[2 * nucl_num + inucl];

				double x = e_coord_0 - n_coord_0;
				double y = e_coord_1 - n_coord_1;
				double z = e_coord_2 - n_coord_2;

				double r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				// Beginning of ao_polynomial computation (now inlined)
				double Y1, Y2, Y3;
				double xy, yz, xz;
				int c, n;
				double da, db, dc, dd;

				// Already computed outsite of the ao_polynomial part
				Y1 = x;
				Y2 = y;
				Y3 = z;

				int llmax = nucleus_max_ang_mom[inucl];
				if (llmax == 0) {
					poly_vgl[0] = 1.;
					poly_vgl[1] = 0.;
					poly_vgl[2] = 0.;
					poly_vgl[3] = 0.;
					poly_vgl[4] = 0.;

					int n = 0;
				} else if (llmax > 0) {
					// Reset pows to 0 for safety. Then we will write over
					// the top left submatrix of size (llmax+3)x3. We will
					// compute indices with llmax and not lmax, so we will
					// use the (llmax+3)*3 first elements of the array
					for (int i = 0; i < 3 * (lmax + 3); i++) {
						pows[i] = 0.;
					}

					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							pows[i + (llmax + 3) * j] = 1.;
						}
					}

					for (int i = 3; i < llmax + 3; i++) {
						pows[i] = pows[(i - 1)] * Y1;
						pows[i + (llmax + 3)] =
							pows[(i - 1) + (llmax + 3)] * Y2;
						pows[i + 2 * (llmax + 3)] =
							pows[(i - 1) + 2 * (llmax + 3)] * Y3;
					}

					for (int i = 0; i < 5; i++) {
						for (int j = 0; j < 4; j++) {
							poly_vgl[i + 5 * j] = 0.;
						}
					}

					poly_vgl[0] = 1.;

					poly_vgl[5] = pows[3];
					poly_vgl[6] = 1.;

					poly_vgl[10] = pows[3 + (llmax + 3)];
					poly_vgl[12] = 1.;

					poly_vgl[15] = pows[3 + 2 * (llmax + 3)];
					poly_vgl[18] = 1.;

					n = 3;
				}

				// l>=2
				dd = 2.;
				for (int d = 2; d <= llmax; d++) {

					da = dd;
					for (int a = d; a >= 0; a--) {

						db = dd - da;

						for (int b = d - a; b >= 0; b--) {

							int c = d - a - b;
							dc = dd - da - db;
							n = n + 1;

							xy = pows[(a + 2)] * pows[(b + 2) + (llmax + 3)];
							yz = pows[(b + 2) + (llmax + 3)] *
								 pows[(c + 2) + 2 * (llmax + 3)];
							xz =
								pows[(a + 2)] * pows[(c + 2) + 2 * (llmax + 3)];

							poly_vgl[5 * (n)] =
								xy * pows[c + 2 + 2 * (llmax + 3)];

							xy = dc * xy;
							xz = db * xz;
							yz = da * yz;

							poly_vgl[1 + 5 * n] = pows[a + 1] * yz;
							poly_vgl[2 + 5 * n] =
								pows[b + 1 + (llmax + 3)] * xz;
							poly_vgl[3 + 5 * n] =
								pows[c + 1 + 2 * (llmax + 3)] * xy;

							poly_vgl[4 + 5 * n] =
								(da - 1.) * pows[a] * yz +
								(db - 1.) * pows[b + (llmax + 3)] * xz +
								(dc - 1.) * pows[c + 2 * (llmax + 3)] * xy;

							db -= 1.;
						}
						da -= 1.;
					}
					dd += 1.;
				}
				// End of ao_polynomial computation (now inlined)
				// poly_vgl is now set from here

				int ishell_start = nucleus_index[inucl];
				int ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;

				// Loop over shells
				int k, l;
				for (int ishell = ishell_start; ishell <= ishell_end;
					 ishell++) {
					k = ao_index[ishell] - 1;
					l = shell_ang_mom[ishell];

					for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2;
						 il++) {

						// value
						ao_vgl[k + 0 * ao_num + ipoint * 5 * ao_num] =
							poly_vgl[il * 5 + 0] *
							shell_vgl[ishell + 0 * shell_num +
									  ipoint * shell_num * 5] *
							ao_factor[k];

						// Grad x
						ao_vgl[k + 1 * ao_num + ipoint * 5 * ao_num] =
							(poly_vgl[il * 5 + 1] *
								 shell_vgl[ishell + 0 * shell_num +
										   ipoint * shell_num * 5] +
							 poly_vgl[il * 5 + 0] *
								 shell_vgl[ishell + 1 * shell_num +
										   ipoint * shell_num * 5]) *
							ao_factor[k];

						// grad y
						ao_vgl[k + 2 * ao_num + ipoint * 5 * ao_num] =
							(poly_vgl[il * 5 + 2] *
								 shell_vgl[ishell + 0 * shell_num +
										   ipoint * shell_num * 5] +
							 poly_vgl[il * 5 + 0] *
								 shell_vgl[ishell + 2 * shell_num +
										   ipoint * shell_num * 5]) *
							ao_factor[k];

						// grad z
						ao_vgl[k + 3 * ao_num + ipoint * 5 * ao_num] =
							(poly_vgl[il * 5 + 3] *
								 shell_vgl[ishell + 0 * shell_num +
										   ipoint * shell_num * 5] +
							 poly_vgl[il * 5 + 0] *
								 shell_vgl[ishell + 3 * shell_num +
										   ipoint * shell_num * 5]) *
							ao_factor[k];

						// Lapl_z
						ao_vgl[k + 4 * ao_num + ipoint * 5 * ao_num] =
							(poly_vgl[il * 5 + 4] *
								 shell_vgl[ishell + 0 * shell_num +
										   ipoint * shell_num * 5] +
							 poly_vgl[il * 5 + 0] *
								 shell_vgl[ishell + 4 * shell_num +
										   ipoint * shell_num * 5] +
							 2.0 * (poly_vgl[il * 5 + 1] *
										shell_vgl[ishell + 1 * shell_num +
												  ipoint * shell_num * 5] +
									poly_vgl[il * 5 + 2] *
										shell_vgl[ishell + 2 * shell_num +
												  ipoint * shell_num * 5] +
									poly_vgl[il * 5 + 3] *
										shell_vgl[ishell + 3 * shell_num +
												  ipoint * shell_num * 5])) *
							ao_factor[k];
						k = k + 1;
					}
				}
			}
		}
		// End of outer compute loop
	}
	// End of target data region
	qmckl_free_device(context, lstart);
	qmckl_free_device(context, poly_vgl_shared);
	qmckl_free_device(context, ao_index);

	qmckl_free_device(context, pows_shared);

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_value */

qmckl_exit_code_device qmckl_compute_ao_value_gaussian_device(
	const qmckl_context_device context, const int64_t ao_num,
	const int64_t shell_num, const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
	const int64_t *restrict nucleus_shell_num, const double *nucleus_range,
	const int32_t *restrict nucleus_max_ang_mom,
	const int32_t *restrict shell_ang_mom, const double *restrict ao_factor,
	double *shell_vgl, double *restrict const ao_value) {

	int64_t n_poly;
	int64_t *lstart;
	double cutoff = 27.631021115928547;

	double *poly_vgl_shared;
	int64_t *powers;
	int64_t *ao_index;

	qmckl_exit_code_device rc;

	qmckl_memory_info_struct_device info;

	lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Multiply "normal" size by point_num to affect subarrays to each thread
	poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * 5 * ao_num * point_num);
	ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);

	// Specific calling function
	int lmax = -1;
#pragma omp target is_device_ptr(nucleus_max_ang_mom) map(tofrom : lmax)
	{
		for (int i = 0; i < nucl_num; i++) {
			if (lmax < nucleus_max_ang_mom[i]) {
				lmax = nucleus_max_ang_mom[i];
			}
		}
	}
	// Multiply "normal" size by point_num to affect subarrays to each thread
	double *pows_shared = qmckl_malloc_device(
		context, sizeof(double) * (lmax + 3) * 3 * point_num);

#pragma omp target is_device_ptr(lstart)
	{
		for (int l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	int k = 1;
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
									 shell_ang_mom, ao_index, lstart)          \
	map(tofrom : k)
	{
		for (int inucl = 0; inucl < nucl_num; inucl++) {
			int ishell_start = nucleus_index[inucl];
			int ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
			for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
				int l = shell_ang_mom[ishell];
				ao_index[ishell] = k;
				k = k + lstart[l + 1] - lstart[l];
			}
		}
	}

#pragma omp target is_device_ptr(                                              \
		ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,     \
			nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared,      \
			nucl_coord, pows_shared, shell_ang_mom, nucleus_range)
	{
		// #pragma omp teams distribute parallel for
		for (int ipoint = 0; ipoint < point_num; ipoint++) {

			// Compute addresses of subarrays from ipoint
			// This way, each thread can write to its own poly_vgl and pows
			// without any race condition
			double *poly_vgl = poly_vgl_shared + ipoint * 5 * ao_num;
			double *pows = pows_shared + ipoint * (lmax + 3) * 3;

			double e_coord_0 = coord[0 * point_num + ipoint];
			double e_coord_1 = coord[1 * point_num + ipoint];
			double e_coord_2 = coord[2 * point_num + ipoint];

			for (int inucl = 0; inucl < nucl_num; inucl++) {

				double n_coord_0 = nucl_coord[0 * nucl_num + inucl];
				double n_coord_1 = nucl_coord[1 * nucl_num + inucl];
				double n_coord_2 = nucl_coord[2 * nucl_num + inucl];

				double x = e_coord_0 - n_coord_0;
				double y = e_coord_1 - n_coord_1;
				double z = e_coord_2 - n_coord_2;

				double r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				// Beginning of ao_polynomial computation (now inlined)
				double Y1, Y2, Y3;
				double xy, yz, xz;
				int c, n;
				double da, db, dc, dd;

				// Already computed outsite of the ao_polynomial part
				Y1 = x;
				Y2 = y;
				Y3 = z;

				int llmax = nucleus_max_ang_mom[inucl];
				if (llmax == 0) {
					poly_vgl[0] = 1.;
					poly_vgl[1] = 0.;
					poly_vgl[2] = 0.;
					poly_vgl[3] = 0.;
					poly_vgl[4] = 0.;

					int n = 0;
				} else if (llmax > 0) {
					// Reset pows to 0 for safety. Then we will write over
					// the top left submatrix of size (llmax+3)x3. We will
					// compute indices with llmax and not lmax, so we will
					// use the (llmax+3)*3 first elements of the array
					for (int i = 0; i < 3 * (lmax + 3); i++) {
						pows[i] = 0.;
					}

					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							pows[i + (llmax + 3) * j] = 1.;
						}
					}

					for (int i = 3; i < llmax + 3; i++) {
						pows[i] = pows[(i - 1)] * Y1;
						pows[i + (llmax + 3)] =
							pows[(i - 1) + (llmax + 3)] * Y2;
						pows[i + 2 * (llmax + 3)] =
							pows[(i - 1) + 2 * (llmax + 3)] * Y3;
					}

					for (int i = 0; i < 5; i++) {
						for (int j = 0; j < 4; j++) {
							poly_vgl[i + 5 * j] = 0.;
						}
					}

					poly_vgl[0] = 1.;

					poly_vgl[5] = pows[3];
					poly_vgl[6] = 1.;

					poly_vgl[10] = pows[3 + (llmax + 3)];
					poly_vgl[12] = 1.;

					poly_vgl[15] = pows[3 + 2 * (llmax + 3)];
					poly_vgl[18] = 1.;

					n = 3;
				}

				// l>=2
				dd = 2.;
				for (int d = 2; d <= llmax; d++) {

					da = dd;
					for (int a = d; a >= 0; a--) {

						db = dd - da;
						for (int b = d - a; b >= 0; b--) {

							int c = d - a - b;
							dc = dd - da - db;
							n = n + 1;

							xy = pows[(a + 2)] * pows[(b + 2) + (llmax + 3)];
							yz = pows[(b + 2) + (llmax + 3)] *
								 pows[(c + 2) + 2 * (llmax + 3)];
							xz =
								pows[(a + 2)] * pows[(c + 2) + 2 * (llmax + 3)];

							poly_vgl[5 * (n)] =
								xy * pows[c + 2 + 2 * (llmax + 3)];

							xy = dc * xy;
							xz = db * xz;
							yz = da * yz;

							poly_vgl[1 + 5 * n] = pows[a + 1] * yz;
							poly_vgl[2 + 5 * n] =
								pows[b + 1 + (llmax + 3)] * xz;
							poly_vgl[3 + 5 * n] =
								pows[c + 1 + 2 * (llmax + 3)] * xy;

							poly_vgl[4 + 5 * n] =
								(da - 1.) * pows[a] * yz +
								(db - 1.) * pows[b + (llmax + 3)] * xz +
								(dc - 1.) * pows[c + 2 * (llmax + 3)] * xy;

							db -= 1.;
						}
						da -= 1.;
					}
					dd += 1.;
				}
				// End of ao_polynomial computation (now inlined)
				// poly_vgl is now set from here

				int ishell_start = nucleus_index[inucl];
				int ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;

				// Loop over shells
				for (int ishell = ishell_start; ishell <= ishell_end;
					 ishell++) {
					int k = ao_index[ishell] - 1;
					int l = shell_ang_mom[ishell];

					for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2;
						 il++) {

						// value
						ao_value[k + ipoint * ao_num] =
							poly_vgl[il * 5 + 0] *
							shell_vgl[ishell + 0 * shell_num +
									  ipoint * shell_num * 5] *
							ao_factor[k];

						k = k + 1;
					}
				}
			}
		}
		// End of outer compute loop
	}
	// End of target data region
	qmckl_free_device(context, lstart);
	qmckl_free_device(context, poly_vgl_shared);
	qmckl_free_device(context, ao_index);

	qmckl_free_device(context, pows_shared);

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// PROVIDE
//**********

/* ao_value */

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_value_device(qmckl_context_device context) {

	qmckl_exit_code_device rc = QMCKL_SUCCESS_DEVICE;

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_ao_basis_ao_value", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_ao_basis_ao_value", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_value_date) {

		/* Allocate array */
		if (ctx->ao_basis.ao_value == NULL) {

			double *ao_value = (double *)qmckl_malloc_device(
				context,
				ctx->ao_basis.ao_num * ctx->point.num * sizeof(double));

			if (ao_value == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_provide_ao_basis_ao_value", NULL);
			}
			ctx->ao_basis.ao_value = ao_value;
		}

		if (ctx->point.date <= ctx->ao_basis.ao_vgl_date &&
			ctx->ao_basis.ao_vgl != NULL) {
			// ao_vgl is already computed and recent enough, we just need to
			// copy the required data to ao_value

			double *v = ctx->ao_basis.ao_value;
			double *vgl = ctx->ao_basis.ao_vgl;
			int point_num = ctx->point.num;
			int ao_num = ctx->ao_basis.ao_num;

#pragma omp target is_device_ptr(v, vgl)
			{
				for (int i = 0; i < point_num; ++i) {
					for (int k = 0; k < ao_num; ++k) {
						v[i * ao_num + k] = vgl[i * ao_num * 5 + k];
					}
				}
			}

		} else {
			// We don't have ao_vgl, so we will compute the values only

			/* Checking for shell_vgl */
			if (ctx->ao_basis.shell_vgl == NULL ||
				ctx->point.date > ctx->ao_basis.shell_vgl_date) {
				qmckl_provide_ao_basis_shell_vgl_device(context);
			}

			if (ctx->ao_basis.type == 'G') {
				rc = qmckl_compute_ao_value_gaussian_device(
					context, ctx->ao_basis.ao_num, ctx->ao_basis.shell_num,
					ctx->point.num, ctx->nucleus.num, ctx->point.coord.data,
					ctx->nucleus.coord.data, ctx->ao_basis.nucleus_index,
					ctx->ao_basis.nucleus_shell_num,
					ctx->ao_basis.nucleus_range,
					ctx->ao_basis.nucleus_max_ang_mom,
					ctx->ao_basis.shell_ang_mom, ctx->ao_basis.ao_factor,
					ctx->ao_basis.shell_vgl, ctx->ao_basis.ao_value);
			} else {
				return qmckl_failwith_device(context, QMCKL_ERRNO_DEVICE,
											 "qmckl_ao_basis_ao_value", NULL);
			}
		}
	}

	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	ctx->ao_basis.ao_value_date = ctx->date;

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// FINALIZE AO BASIS
//**********

qmckl_exit_code_device
qmckl_finalize_ao_basis_hpc_device(qmckl_context_device context) {

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;

	int device_id = qmckl_get_device_id(context);

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

#pragma omp target map(tofrom : shell_max_ptr[ : 1], prim_max_ptr[ : 1])       \
	is_device_ptr(nucleus_shell_num, nucleus_index, shell_prim_num,            \
					  prim_num_per_nucleus)
	{

		for (int inucl = 0; inucl < nucl_num; ++inucl) {
			shell_max_ptr[0] = nucleus_shell_num[inucl] > shell_max_ptr[0]
								   ? nucleus_shell_num[inucl]
								   : shell_max_ptr[0];

			int64_t prim_num = 0;
			int64_t ishell_start = nucleus_index[inucl];
			int64_t ishell_end =
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

#pragma omp target is_device_ptr(                                              \
		expo_expo, expo_index, coef, newcoef, nucleus_index, shell_prim_index, \
			nucleus_shell_num, exponent, coefficient_normalized,               \
			shell_prim_num, expo_per_nucleus_data, coef_per_nucleus_data,      \
			prim_num_per_nucleus, newidx)
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
			int64_t ishell_start = nucleus_index[inucl];
			int64_t ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl];

			for (int64_t ishell = ishell_start; ishell < ishell_end; ++ishell) {

				int64_t iprim_start = shell_prim_index[ishell];
				int64_t iprim_end =
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
				int64_t iprim_start = shell_prim_index[ishell];
				int64_t iprim_end =
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

	qmckl_free_device(context, expo_expo);
	qmckl_free_device(context, expo_index);
	qmckl_free_device(context, coef);
	qmckl_free_device(context, newcoef);
	qmckl_free_device(context, newidx);

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_finalize_ao_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_finalize_ao_basis_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

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

#pragma omp target is_device_ptr(shell_prim_index, shell_prim_num,             \
									 coefficient_normalized, coefficient,      \
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

#pragma omp target is_device_ptr(nucleus_range, nucleus_index,                 \
									 nucleus_shell_num, shell_prim_index,      \
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
