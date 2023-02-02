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
	cutoff = 27.631021115928547; //-dlog(1.d-12)

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

qmckl_exit_code qmckl_compute_ao_vgl_gaussian_device(
	const qmckl_context context, const int64_t ao_num, const int64_t shell_num,
	const int64_t point_num, const int64_t nucl_num,
	const double *restrict coord, const double *restrict nucl_coord,
	const int64_t *restrict nucleus_index,
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
	double cutoff = 27.631021115928547;
	int64_t qmckl_ao_polynomial_vgl_doc_f;
	int64_t size_max = 0;

	double *poly_vgl;
	int64_t *powers;
	int64_t *ao_index;

	qmckl_exit_code rc;
	int lmax, c;

	qmckl_memory_info_struct info;

	info.size = sizeof(int64_t) * 21;
	lstart = qmckl_malloc_device(context, info);

	info.size = sizeof(double) * 8;
	e_coord = qmckl_malloc_device(context, info);
	info.size = sizeof(double) * 8;
	n_coord = qmckl_malloc_device(context, info);

	info.size = sizeof(double) * 5 * ao_num;
	poly_vgl = qmckl_malloc_device(context, info);
	info.size = sizeof(int64_t) * ao_num;
	ao_index = qmckl_malloc_device(context, info);

	// Specific calling function
	lmax = -1;
#pragma omp target is_device_ptr(nucleus_max_ang_mom)
	{
#pragma omp target update to(lmax)
		{
			for (int i = 0; i < nucl_num; i++) {
				if (lmax < nucleus_max_ang_mom[i]) {
					lmax = nucleus_max_ang_mom[i];
				}
			}
		}
	}
	info.size = sizeof(double) * (lmax + 3) * 3;
	double *pows = qmckl_malloc_device(context, info);

#pragma omp target is_device_ptr(lstart)
	{
		for (l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	k = 1;
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
								 shell_ang_mom, ao_index, lstart)
	{
#pragma omp target update to(k)
		{
			for (inucl = 0; inucl < nucl_num; inucl++) {
				ishell_start = nucleus_index[inucl];
				ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
				for (ishell = ishell_start; ishell <= ishell_end; ishell++) {
					l = shell_ang_mom[ishell];
					ao_index[ishell] = k;
					k = k + lstart[l + 1] - lstart[l];
				}
			}
		}
	}

#pragma omp target is_device_ptr(                                              \
	ao_vgl, lstart, e_coord, n_coord, ao_index, ao_factor, coord,              \
	nucleus_max_ang_mom, nucleus_index, nucleus_shell_num, shell_vgl,          \
	poly_vgl, nucl_coord, pows, shell_ang_mom, nucleus_range)
	{
		//#pragma omp teams distribute parallel for
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

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				// Beginning of ao_polynomial computation (now inlined)
				double Y1, Y2, Y3;
				double xy, yz, xz;
				int c, n;
				double da, db, dc, dd;

				Y1 = e_coord[0] - n_coord[0];
				Y2 = e_coord[1] - n_coord[1];
				Y3 = e_coord[2] - n_coord[2];

				int llmax = nucleus_max_ang_mom[inucl];
				if (llmax == 0) {
					poly_vgl[0] = 1.;
					poly_vgl[1] = 0.;
					poly_vgl[2] = 0.;
					poly_vgl[3] = 0.;
					poly_vgl[4] = 0.;

					int n = 0;
				} else if (llmax > 0) {
					// Reset pows to 0 for safety. Then we will write over the
					// top left submatrix of size (llmax+3)x3. We will compute
					// indices with llmax and not lmax, so we will use the
					// (llmax+3)*3 first elements of the array
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

				ishell_start = nucleus_index[inucl];
				ishell_end =
					nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;

				// Loop over shells
				for (ishell = ishell_start; ishell <= ishell_end; ishell++) {
					k = ao_index[ishell] - 1;
					l = shell_ang_mom[ishell];

					for (il = lstart[l] - 1; il <= lstart[l + 1] - 2; il++) {

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
	}
	qmckl_free_device(context, lstart);
	qmckl_free_device(context, e_coord);
	qmckl_free_device(context, n_coord);
	qmckl_free_device(context, poly_vgl);
	qmckl_free_device(context, ao_index);

	qmckl_free_device(context, pows);

	return QMCKL_SUCCESS;
}
