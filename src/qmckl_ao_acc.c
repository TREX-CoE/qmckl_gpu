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

	int iprim_start, iprim_end;
	double x, y, z, two_a, ar2, r2, v, cutoff;

	qmckl_exit_code info = QMCKL_SUCCESS;

	// Don't compute exponentials when the result will be almost zero.
	// TODO : Use numerical precision here
	cutoff = 27.631021115928547; //-dlog(1.d-12)

    int* shell_to_nucl = qmckl_malloc_device(context, sizeof(int)*shell_num);

#pragma acc data deviceptr(nucleus_shell_num, nucleus_index, nucleus_range,    \
							   shell_prim_index, shell_prim_num, coord,        \
							   nucl_coord, expo, coef_normalized, shell_vgl,    \
                               shell_to_nucl)
	{

#pragma acc kernels
    { 
    for (int inucl = 0; inucl < nucl_num; inucl++) {
        int ishell_start = nucleus_index[inucl];
        int ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
        for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
            shell_to_nucl[ishell] = inucl ;
        }
    }
    }
		/*
		 * BUG As of now, "acc parallel" caused the following internal compiler
		 * error on  gcc (Spack GCC) 12.1.0 :
		 *
		 * ../src/qmckl_ao_acc.c: In function
		 * 'qmckl_compute_ao_basis_shell_gaussian_vgl_device._omp_fn.0':
		 * ../src/qmckl_ao_acc.c:31:9: internal compiler error: in
		 * expand_UNIQUE, at internal-fn.cc:2996 31 | #pragma acc parallel loop
		 * gang worker vector
		 *
		 *  TODO Until this error is fixed, we might want to wrap desired
		 * pragmas in #ifdefs depending on the compiler 
		 * */
#pragma acc parallel loop collapse(2) 
		for (int ipoint = 0; ipoint < point_num; ipoint++) {
            for (int ishell = 0; ishell < shell_num; ishell++) {

			    int inucl = shell_to_nucl[ishell];

				x = coord[ipoint] - nucl_coord[inucl];
				y = coord[ipoint + point_num] - nucl_coord[inucl + nucl_num];
				z = coord[ipoint + 2 * point_num] -
					nucl_coord[inucl + 2 * nucl_num];

				r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				double t0 = 0;
				double t1 = 0;
				double t2 = 0;
				double t3 = 0;
				double t4 = 0;

				iprim_start = shell_prim_index[ishell];
				iprim_end = shell_prim_index[ishell] + shell_prim_num[ishell] - 1;

#pragma acc loop seq
				for (int iprim = iprim_start; iprim <= iprim_end; iprim++) {

					ar2 = expo[iprim] * r2;
					if (ar2 > cutoff) {
						continue;
					}

					v = coef_normalized[iprim] * exp(-ar2);
					two_a = -2 * expo[iprim] * v;

					t0+= v; 
					t1+= two_a * x;
					t2+= two_a * y;
					t3+= two_a * z;
					t4+= two_a * (3 - 2 * ar2);
                    
				}

				shell_vgl[ishell + 0 * shell_num +
						  ipoint * shell_num * 5] = t0;

				shell_vgl[ishell + 1 * shell_num +
						  ipoint * shell_num * 5] = t1;

				shell_vgl[ishell + 2 * shell_num +
						  ipoint * shell_num * 5] = t2;

				shell_vgl[ishell + 3 * shell_num +
						  ipoint * shell_num * 5] = t3;

				shell_vgl[ishell + 4 * shell_num +
					  ipoint * shell_num * 5] = t4; 
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

	int64_t n_poly;
	int64_t *lstart;
	double cutoff = 27.631021115928547;

	double *poly_vgl_shared;
	int64_t *powers;
	int64_t *ao_index;

	qmckl_exit_code rc;

	lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Multiply "normal" size by point_num to affect subarrays to each thread
	poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * 5 * ao_num * point_num);
	ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);

	// Specific calling function
	int lmax = -1;
	int *lmax_p = &lmax;
#pragma acc data deviceptr(nucleus_max_ang_mom) copyin(lmax_p[0 : 1])
	{
#pragma acc kernels
		{
			for (int i = 0; i < nucl_num; i++) {
				if (lmax_p[0] < nucleus_max_ang_mom[i]) {
					lmax_p[0] = nucleus_max_ang_mom[i];
				}
			}
		}
#pragma acc update host(lmax_p[0 : 1])
	}
	// Multiply "normal" size by point_num to affect subarrays to each thread
	double *pows_shared = qmckl_malloc_device(
		context, sizeof(double) * (lmax + 3) * 3 * point_num);

#pragma acc kernels deviceptr(lstart)
	{
		for (int l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	int k = 1;
	int *k_p = &k;
#pragma acc data deviceptr(nucleus_index, nucleus_shell_num, shell_ang_mom,    \
							   ao_index, lstart) copyin(k_p[0 : 1])
	{
#pragma acc kernels
		{for (int inucl = 0; inucl < nucl_num;
			  inucl++){int ishell_start = nucleus_index[inucl];
	int ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
	for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
		int l = shell_ang_mom[ishell];
		ao_index[ishell] = k_p[0];
		k_p[0] = k_p[0] + lstart[l + 1] - lstart[l];
	}
}
}
#pragma acc update host(k_p[0 : 1])
}

#pragma acc data deviceptr(                                                    \
		ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,       \
			nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared,      \
			nucl_coord, pows_shared, shell_ang_mom, nucleus_range)
{

    double poly_vgl[5*ao_num][point_num];

	// BUG See qmckl_compute_ao_basis_shell_gaussian_vgl_device above
#pragma acc parallel loop gang worker vector
	for (int ipoint = 0; ipoint < point_num; ipoint++) {

		// Compute addresses of subarrays from ipoint
		// This way, each thread can write to its own poly_vgl and pows
		// without any race condition
		//double *poly_vgl = poly_vgl_shared + ipoint * 5 * ao_num;
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
				poly_vgl[0][ipoint] = 1.;
				poly_vgl[1][ipoint] = 0.;
				poly_vgl[2][ipoint] = 0.;
				poly_vgl[3][ipoint] = 0.;
				poly_vgl[4][ipoint] = 0.;

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
					pows[i + (llmax + 3)] = pows[(i - 1) + (llmax + 3)] * Y2;
					pows[i + 2 * (llmax + 3)] =
						pows[(i - 1) + 2 * (llmax + 3)] * Y3;
				}

				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 4; j++) {
						poly_vgl[i + 5 * j][ipoint] = 0.;
					}
				}

				poly_vgl[0][ipoint] = 1.;

				poly_vgl[5][ipoint] = pows[3];
				poly_vgl[6][ipoint] = 1.;

				poly_vgl[10][ipoint] = pows[3 + (llmax + 3)];
				poly_vgl[12][ipoint] = 1.;

				poly_vgl[15][ipoint] = pows[3 + 2 * (llmax + 3)];
				poly_vgl[18][ipoint] = 1.;

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
						xz = pows[(a + 2)] * pows[(c + 2) + 2 * (llmax + 3)];

						poly_vgl[5 * (n)][ipoint] = xy * pows[c + 2 + 2 * (llmax + 3)];

						xy = dc * xy;
						xz = db * xz;
						yz = da * yz;

						poly_vgl[1 + 5 * n][ipoint] = pows[a + 1] * yz;
						poly_vgl[2 + 5 * n][ipoint] = pows[b + 1 + (llmax + 3)] * xz;
						poly_vgl[3 + 5 * n][ipoint] =
							pows[c + 1 + 2 * (llmax + 3)] * xy;

						poly_vgl[4 + 5 * n][ipoint] =
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
			for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
				int k = ao_index[ishell] - 1;
				int l = shell_ang_mom[ishell];

				for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2; il++) {

					// value
					ao_vgl[k + 0 * ao_num + ipoint * 5 * ao_num] =
						poly_vgl[il * 5 + 0][ipoint] *
						shell_vgl[ishell + 0 * shell_num +
								  ipoint * shell_num * 5] *
						ao_factor[k];

					// Grad x
					ao_vgl[k + 1 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 1][ipoint] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][ipoint] *
							 shell_vgl[ishell + 1 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// grad y
					ao_vgl[k + 2 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 2][ipoint] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][ipoint] *
							 shell_vgl[ishell + 2 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// grad z
					ao_vgl[k + 3 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 3][ipoint] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][ipoint] *
							 shell_vgl[ishell + 3 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// Lapl_z
					ao_vgl[k + 4 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 4][ipoint] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][ipoint] *
							 shell_vgl[ishell + 4 * shell_num +
									   ipoint * shell_num * 5] +
						 2.0 * (poly_vgl[il * 5 + 1][ipoint] *
									shell_vgl[ishell + 1 * shell_num +
											  ipoint * shell_num * 5] +
								poly_vgl[il * 5 + 2][ipoint] *
									shell_vgl[ishell + 2 * shell_num +
											  ipoint * shell_num * 5] +
								poly_vgl[il * 5 + 3][ipoint] *
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

return QMCKL_SUCCESS;
}

/* ao_value */

qmckl_exit_code qmckl_compute_ao_value_gaussian_device(
	const qmckl_context context, const int64_t ao_num, const int64_t shell_num,
	const int64_t point_num, const int64_t nucl_num,
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

	qmckl_exit_code rc;

	lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Multiply "normal" size by point_num to affect subarrays to each thread
	poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * 5 * ao_num * point_num);
	ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);

	// Specific calling function
	int lmax = -1;
	int *lmax_p = &lmax;
#pragma acc data deviceptr(nucleus_max_ang_mom) copyin(lmax_p[0 : 1])
	{
#pragma acc kernels
		{
			for (int i = 0; i < nucl_num; i++) {
				if (lmax_p[0] < nucleus_max_ang_mom[i]) {
					lmax_p[0] = nucleus_max_ang_mom[i];
				}
			}
		}
#pragma acc update host(lmax_p[0 : 1])
	}
	// Multiply "normal" size by point_num to affect subarrays to each thread
	double *pows_shared = qmckl_malloc_device(
		context, sizeof(double) * (lmax + 3) * 3 * point_num);

#pragma acc kernels deviceptr(lstart)
	{
		for (int l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}

	int k = 1;
	int *k_p = &k;
#pragma acc data deviceptr(nucleus_index, nucleus_shell_num, shell_ang_mom,    \
							   ao_index, lstart) copyin(k_p[0 : 1])
	{
#pragma acc kernels
		{for (int inucl = 0; inucl < nucl_num;
			  inucl++){int ishell_start = nucleus_index[inucl];
	int ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
	for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
		int l = shell_ang_mom[ishell];
		ao_index[ishell] = k_p[0];
		k_p[0] = k_p[0] + lstart[l + 1] - lstart[l];
	}
}
}
#pragma acc update host(k_p[0 : 1])
}

#pragma acc data deviceptr(                                                    \
		ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,     \
			nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared,      \
			nucl_coord, pows_shared, shell_ang_mom, nucleus_range)
{

	// BUG See qmckl_compute_ao_basis_shell_gaussian_vgl_device above

#pragma acc parallel loop gang worker vector
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
					pows[i + (llmax + 3)] = pows[(i - 1) + (llmax + 3)] * Y2;
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
						xz = pows[(a + 2)] * pows[(c + 2) + 2 * (llmax + 3)];

						poly_vgl[5 * (n)] = xy * pows[c + 2 + 2 * (llmax + 3)];

						xy = dc * xy;
						xz = db * xz;
						yz = da * yz;

						poly_vgl[1 + 5 * n] = pows[a + 1] * yz;
						poly_vgl[2 + 5 * n] = pows[b + 1 + (llmax + 3)] * xz;
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
			for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
				int k = ao_index[ishell] - 1;
				int l = shell_ang_mom[ishell];

				for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2; il++) {

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

return QMCKL_SUCCESS;
}

//**********
// PROVIDE
//**********

/* ao_value */

qmckl_exit_code qmckl_provide_ao_basis_ao_value_device(qmckl_context context) {

	qmckl_exit_code rc = QMCKL_SUCCESS;

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
		return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
							  "qmckl_provide_ao_basis_ao_value", NULL);
	}

	qmckl_context_struct *const ctx = (qmckl_context_struct *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith(context, QMCKL_NOT_PROVIDED,
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
				return qmckl_failwith(context, QMCKL_ALLOCATION_FAILED,
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

#pragma acc kernels deviceptr(v, vgl)
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
				return qmckl_failwith(context, QMCKL_ERRNO,
									  "qmckl_ao_basis_ao_value", NULL);
			}
		}
	}

	if (rc != QMCKL_SUCCESS) {
		return rc;
	}

	ctx->ao_basis.ao_value_date = ctx->date;

	return QMCKL_SUCCESS;
}
