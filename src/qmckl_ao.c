#include "include/qmckl_ao.h"
#define MAX_MEMORY_SIZE 16.0 * 1024 * 1024 * 1024

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

	// Don't compute exponentials when the result will be almost zero.
	// TODO : Use numerical precision here
	double cutoff = 27.631021115928547; //-dlog(1.d-12)

	int *shell_to_nucl = qmckl_malloc_device(context, sizeof(int) * shell_num);

#pragma acc parallel loop \
            deviceptr(shell_to_nucl, nucleus_index, nucleus_shell_num)
#pragma omp target teams distribute parallel for \
            is_device_ptr(shell_to_nucl, nucleus_index, nucleus_shell_num)
		for (int inucl = 0; inucl < nucl_num; inucl++) {
			int ishell_start = nucleus_index[inucl];
			int ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
			for (int ishell = ishell_start; ishell <= ishell_end;
				 ishell++) {
				shell_to_nucl[ishell] = inucl;
			}
		}
	
#pragma acc parallel loop collapse(2) \
            deviceptr(nucleus_shell_num, nucleus_index, nucleus_range,    \
            		  shell_prim_index, shell_prim_num, coord,            \
            		  nucl_coord, expo, coef_normalized, shell_vgl,       \
            		  shell_to_nucl)
#pragma omp target teams loop collapse(2) \
            is_device_ptr(nucleus_shell_num, nucleus_index, nucleus_range,    \
            		      shell_prim_index, shell_prim_num, coord,            \
            		      nucl_coord, expo, coef_normalized, shell_vgl,       \
            		      shell_to_nucl)
		for (int ipoint = 0; ipoint < point_num; ipoint++) {
			for (int ishell = 0; ishell < shell_num; ishell++) {

				int inucl = shell_to_nucl[ishell];

				double x = coord[ipoint] - nucl_coord[inucl];
				double y =
					coord[ipoint + point_num] - nucl_coord[inucl + nucl_num];
				double z = coord[ipoint + 2 * point_num] -
						   nucl_coord[inucl + 2 * nucl_num];

				double r2 = x * x + y * y + z * z;

				if (r2 > cutoff * nucleus_range[inucl]) {
					continue;
				}

				double t0 = 0;
				double t1 = 0;
				double t2 = 0;
				double t3 = 0;
				double t4 = 0;

				int iprim_start = shell_prim_index[ishell];
				int iprim_end =
					shell_prim_index[ishell] + shell_prim_num[ishell] - 1;

				// BEWARE: with omp set noautopar with nvidia compilers,
				// this loop must be executed serially not to degrade performances.
#pragma acc loop seq
				for (int iprim = iprim_start; iprim <= iprim_end; iprim++) {

					double ar2 = expo[iprim] * r2;
					if (ar2 > cutoff) {
						continue;
					}

					double v = coef_normalized[iprim] * exp(-ar2);
					double two_a = -2 * expo[iprim] * v;

					t0 += v;
					t1 += two_a * x;
					t2 += two_a * y;
					t3 += two_a * z;
					t4 += two_a * (3 - 2 * ar2);
				}

				shell_vgl[ishell + 0 * shell_num + ipoint * shell_num * 5] = t0;

				shell_vgl[ishell + 1 * shell_num + ipoint * shell_num * 5] = t1;

				shell_vgl[ishell + 2 * shell_num + ipoint * shell_num * 5] = t2;

				shell_vgl[ishell + 3 * shell_num + ipoint * shell_num * 5] = t3;

				shell_vgl[ishell + 4 * shell_num + ipoint * shell_num * 5] = t4;
			}
		}

	qmckl_free_device(context, shell_to_nucl);

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;
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

	double cutoff = 27.631021115928547;

	// TG: MAX_MEMORY_SIZE should be roughly the GPU RAM, to be set at configure
	// time? Not to exceed GPU memory when allocating poly_vgl
	// int64_t target_chunk = (MAX_MEMORY_SIZE / 4.) / (sizeof(double) * 5 * ao_num);
	// A good guess, can result in large memory occupation though
	int64_t target_chunk = 128*1024;
	size_t max_chunk_size = ((target_chunk) / nucl_num) * nucl_num;
	int64_t num_iters = point_num * nucl_num;
	int64_t chunk_size =
		(num_iters < max_chunk_size) ? num_iters : max_chunk_size;
	int64_t num_sub_iters = (num_iters + chunk_size - 1) / chunk_size;
	int64_t poly_dim = 5 * ao_num * chunk_size;

	double *poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * poly_dim);
	int64_t *ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);
	int64_t *lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Specific calling function
	int lmax = -1;
#pragma acc kernels deviceptr(nucleus_max_ang_mom) copy(lmax)
#pragma omp target is_device_ptr(nucleus_max_ang_mom) map(tofrom : lmax)
	{
		for (int i = 0; i < nucl_num; i++){
			if (lmax < nucleus_max_ang_mom[i]){
				lmax = nucleus_max_ang_mom[i];
			}
		}
	}
	size_t pows_dim = (lmax + 3) * 3 * chunk_size;
	double *pows_shared = qmckl_malloc_device(context, sizeof(double) * pows_dim);

#pragma acc kernels deviceptr(lstart)
#pragma omp target is_device_ptr(lstart)
	{
		for (int l = 0; l < 21; l++) {
			lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
		}
	}
	
	int k = 1;
	int *shell_to_nucl = qmckl_malloc_device(context, sizeof(int) * shell_num);
#pragma acc kernels deviceptr(nucleus_index, nucleus_shell_num, shell_ang_mom,    \
						   ao_index, lstart, shell_to_nucl) copy(k)
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
								 shell_ang_mom, ao_index, lstart,              \
								 shell_to_nucl) map(tofrom : k)
	{
	for (int inucl = 0; inucl < nucl_num; inucl++){
		int ishell_start = nucleus_index[inucl];
		int ishell_end = nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
		for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
			int l = shell_ang_mom[ishell];
			ao_index[ishell] = k;
			k = k + lstart[l + 1] - lstart[l];
			shell_to_nucl[ishell] = inucl;
		}
	}
	}

	double(*poly_vgl)[chunk_size] = (double(*)[chunk_size])poly_vgl_shared;
	double(*pows)[chunk_size] = (double(*)[chunk_size])pows_shared;
	for (int sub_iter = 0; sub_iter < num_sub_iters; sub_iter++) {
#pragma acc parallel loop gang worker vector async(1) \
    deviceptr( \
	ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,           \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl, nucl_coord,  \
	pows_shared, pows, shell_ang_mom, nucleus_range, shell_to_nucl)
#pragma omp target teams distribute parallel for \
    is_device_ptr( \
	ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,           \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl, nucl_coord,  \
	pows_shared, pows, shell_ang_mom, nucleus_range, shell_to_nucl)
		for (int iter = 0; iter < chunk_size; iter++) {

			int step = iter + sub_iter * chunk_size;
			if (step >= num_iters)
				continue;

			int ipoint = step / nucl_num;
			int inucl = step % nucl_num;

			double e_coord_0 = coord[0 * point_num + ipoint];
			double e_coord_1 = coord[1 * point_num + ipoint];
			double e_coord_2 = coord[2 * point_num + ipoint];

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
			int n;

			// Already computed outsite of the ao_polynomial part
			double Y1 = x;
			double Y2 = y;
			double Y3 = z;

			int llmax = nucleus_max_ang_mom[inucl];
			if (llmax == 0) {
				poly_vgl[0][iter] = 1.;
				poly_vgl[1][iter] = 0.;
				poly_vgl[2][iter] = 0.;
				poly_vgl[3][iter] = 0.;
				poly_vgl[4][iter] = 0.;

				n = 0;
			} else if (llmax > 0) {
				// Reset pows to 0 for safety. Then we will write over the
				// top left submatrix of size (llmax+3)x3. We will compute
				// indices with llmax and not lmax, so we will use the
				// (llmax+3)*3 first elements of the array
				for (int i = 0; i < 3 * (lmax + 3); i++) {
					pows[i][iter] = 0.;
				}

				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						pows[i + (llmax + 3) * j][iter] = 1.;
					}
				}

				for (int i = 3; i < llmax + 3; i++) {
					pows[i][iter] = pows[(i - 1)][iter] * Y1;
					pows[i + (llmax + 3)][iter] =
						pows[(i - 1) + (llmax + 3)][iter] * Y2;
					pows[i + 2 * (llmax + 3)][iter] =
						pows[(i - 1) + 2 * (llmax + 3)][iter] * Y3;
				}

				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 4; j++) {
						poly_vgl[i + 5 * j][iter] = 0.;
					}
				}

				poly_vgl[0][iter] = 1.;

				poly_vgl[5][iter] = pows[3][iter];
				poly_vgl[6][iter] = 1.;

				poly_vgl[10][iter] = pows[3 + (llmax + 3)][iter];
				poly_vgl[12][iter] = 1.;

				poly_vgl[15][iter] = pows[3 + 2 * (llmax + 3)][iter];
				poly_vgl[18][iter] = 1.;

				n = 3;
			}

			double xy, yz, xz;
			double da, db, dc, dd;

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

						xy = pows[(a + 2)][iter] *
							 pows[(b + 2) + (llmax + 3)][iter];
						yz = pows[(b + 2) + (llmax + 3)][iter] *
							 pows[(c + 2) + 2 * (llmax + 3)][iter];
						xz = pows[(a + 2)][iter] *
							 pows[(c + 2) + 2 * (llmax + 3)][iter];

						poly_vgl[5 * (n)][iter] =
							xy * pows[c + 2 + 2 * (llmax + 3)][iter];

						xy = dc * xy;
						xz = db * xz;
						yz = da * yz;

						poly_vgl[1 + 5 * n][iter] = pows[a + 1][iter] * yz;
						poly_vgl[2 + 5 * n][iter] =
							pows[b + 1 + (llmax + 3)][iter] * xz;
						poly_vgl[3 + 5 * n][iter] =
							pows[c + 1 + 2 * (llmax + 3)][iter] * xy;

						poly_vgl[4 + 5 * n][iter] =
							(da - 1.) * pows[a][iter] * yz +
							(db - 1.) * pows[b + (llmax + 3)][iter] * xz +
							(dc - 1.) * pows[c + 2 * (llmax + 3)][iter] * xy;

						db -= 1.;
					}
					da -= 1.;
				}
				dd += 1.;
			}
			// End of ao_polynomial computation (now inlined)
			// poly_vgl is now set from here
		}

#pragma acc parallel loop collapse(2) async(1) \
    deviceptr( \
	ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,           \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl, nucl_coord,  \
	pows_shared, pows, shell_ang_mom, nucleus_range, shell_to_nucl)
#pragma omp target teams loop collapse(2) \
    is_device_ptr( \
	ao_vgl, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,           \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl, nucl_coord,  \
	pows_shared, pows, shell_ang_mom, nucleus_range, shell_to_nucl)
		for (int iter_new = 0; iter_new < chunk_size / nucl_num; iter_new++) {
			for (int ishell = 0; ishell < shell_num; ishell++) {

				int ipoint = iter_new + chunk_size / nucl_num * sub_iter;
				int inucl = shell_to_nucl[ishell];
				int iter = iter_new * nucl_num + inucl;

				if (ipoint >= point_num)
					continue;

				double e_coord_0 = coord[0 * point_num + ipoint];
				double e_coord_1 = coord[1 * point_num + ipoint];
				double e_coord_2 = coord[2 * point_num + ipoint];

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

				int k = ao_index[ishell] - 1;
				int l = shell_ang_mom[ishell];

				for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2; il++) {

					// value
					ao_vgl[k + 0 * ao_num + ipoint * 5 * ao_num] =
						poly_vgl[il * 5 + 0][iter] *
						shell_vgl[ishell + 0 * shell_num +
								  ipoint * shell_num * 5] *
						ao_factor[k];

					// Grad x
					ao_vgl[k + 1 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 1][iter] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][iter] *
							 shell_vgl[ishell + 1 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// grad y
					ao_vgl[k + 2 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 2][iter] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][iter] *
							 shell_vgl[ishell + 2 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// grad z
					ao_vgl[k + 3 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 3][iter] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][iter] *
							 shell_vgl[ishell + 3 * shell_num +
									   ipoint * shell_num * 5]) *
						ao_factor[k];

					// Lapl_z (_z?)
					ao_vgl[k + 4 * ao_num + ipoint * 5 * ao_num] =
						(poly_vgl[il * 5 + 4][iter] *
							 shell_vgl[ishell + 0 * shell_num +
									   ipoint * shell_num * 5] +
						 poly_vgl[il * 5 + 0][iter] *
							 shell_vgl[ishell + 4 * shell_num +
									   ipoint * shell_num * 5] +
						 2.0 * (poly_vgl[il * 5 + 1][iter] *
									shell_vgl[ishell + 1 * shell_num +
											  ipoint * shell_num * 5] +
								poly_vgl[il * 5 + 2][iter] *
									shell_vgl[ishell + 2 * shell_num +
											  ipoint * shell_num * 5] +
								poly_vgl[il * 5 + 3][iter] *
									shell_vgl[ishell + 3 * shell_num +
											  ipoint * shell_num * 5])) *
						ao_factor[k];
					k = k + 1;
				}
			}
		}
#pragma acc wait(1)
	// End of outer compute loop
}
// End of target data region
qmckl_free_device(context, ao_index);
qmckl_free_device(context, poly_vgl_shared);
qmckl_free_device(context, pows_shared);
qmckl_free_device(context, shell_to_nucl);
qmckl_free_device(context, lstart);

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

	double cutoff = 27.631021115928547;

	int64_t target_chunk = 128*1024;
	// int64_t target_chunk = (MAX_MEMORY_SIZE / 4.) / (sizeof(double) * ao_num);
	size_t max_chunk_size = ((target_chunk) / nucl_num) * nucl_num;
	int64_t num_iters = point_num * nucl_num;
	int64_t chunk_size =
		(num_iters < max_chunk_size) ? num_iters : max_chunk_size;
	int64_t num_sub_iters = (num_iters + chunk_size - 1) / chunk_size;
	int64_t poly_dim = ao_num * chunk_size;

	double *poly_vgl_shared =
		qmckl_malloc_device(context, sizeof(double) * poly_dim);
	int64_t *ao_index = qmckl_malloc_device(context, sizeof(int64_t) * ao_num);
	int64_t *lstart = qmckl_malloc_device(context, sizeof(int64_t) * 21);

	// Specific calling function
	int lmax = -1;
#pragma acc kernels deviceptr(nucleus_max_ang_mom) copy(lmax)
#pragma omp target is_device_ptr(nucleus_max_ang_mom) map(tofrom : lmax)
	{
		for (int i = 0; i < nucl_num; i++) {
			if (lmax < nucleus_max_ang_mom[i]) {
				lmax = nucleus_max_ang_mom[i];
			}
		}
	}
size_t pows_dim = (lmax + 3) * 3 * chunk_size;
double *pows_shared = qmckl_malloc_device(context, sizeof(double) * pows_dim);

#pragma acc kernels deviceptr(lstart)
#pragma omp target is_device_ptr(lstart)
{
	for (int l = 0; l < 21; l++) {
		lstart[l] = l * (l + 1) * (l + 2) / 6 + 1;
	}
}

int k = 1;
int *shell_to_nucl = qmckl_malloc_device(context, sizeof(int) * shell_num);
#pragma acc kernels deviceptr(nucleus_index, nucleus_shell_num, shell_ang_mom,    \
						   ao_index, lstart, shell_to_nucl) copy(k)
#pragma omp target is_device_ptr(nucleus_index, nucleus_shell_num,             \
								 shell_ang_mom, ao_index, lstart,              \
								 shell_to_nucl) map(tofrom : k)
{
		for (int inucl = 0; inucl < nucl_num; inucl++) {
			int ishell_start = nucleus_index[inucl];
			int ishell_end =
				nucleus_index[inucl] + nucleus_shell_num[inucl] - 1;
			for (int ishell = ishell_start; ishell <= ishell_end; ishell++) {
				int l = shell_ang_mom[ishell];
				ao_index[ishell] = k;
				k = k + lstart[l + 1] - lstart[l];
				shell_to_nucl[ishell] = inucl;
			}
		}

}

	double(*poly_vgl)[chunk_size] = (double(*)[chunk_size])poly_vgl_shared;
	double(*pows)[chunk_size] = (double(*)[chunk_size])pows_shared;
	for (int sub_iter = 0; sub_iter < num_sub_iters; sub_iter++) {
#pragma acc parallel loop gang worker vector async(1) \
	deviceptr( \
	ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,         \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl,    \
	pows, nucl_coord, pows_shared, shell_ang_mom, nucleus_range,               \
	shell_to_nucl)
#pragma omp target teams distribute parallel for \
	is_device_ptr( \
	ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,         \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl,    \
	pows, nucl_coord, pows_shared, shell_ang_mom, nucleus_range,               \
	shell_to_nucl)
		for (int iter = 0; iter < chunk_size; iter++) {

			int step = iter + sub_iter * chunk_size;
			if (step >= num_iters)
				continue;

			int ipoint = step / nucl_num;
			int inucl = step % nucl_num;

			double e_coord_0 = coord[0 * point_num + ipoint];
			double e_coord_1 = coord[1 * point_num + ipoint];
			double e_coord_2 = coord[2 * point_num + ipoint];

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
			int n;

			// Already computed outsite of the ao_polynomial part
			double Y1 = x;
			double Y2 = y;
			double Y3 = z;

			int llmax = nucleus_max_ang_mom[inucl];
			if (llmax == 0) {
				poly_vgl[0][iter] = 1.;
				poly_vgl[1][iter] = 0.;
				poly_vgl[2][iter] = 0.;
				poly_vgl[3][iter] = 0.;
				poly_vgl[4][iter] = 0.;

				n = 0;
			} else if (llmax > 0) {
				// Reset pows to 0 for safety. Then we will write over
				// the top left submatrix of size (llmax+3)x3. We will
				// compute indices with llmax and not lmax, so we will
				// use the (llmax+3)*3 first elements of the array
				for (int i = 0; i < 3 * (lmax + 3); i++) {
					pows[i][iter] = 0.;
				}

				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						pows[i + (llmax + 3) * j][iter] = 1.;
					}
				}

				for (int i = 3; i < llmax + 3; i++) {
					pows[i][iter] = pows[(i - 1)][iter] * Y1;
					pows[i + (llmax + 3)][iter] =
						pows[(i - 1) + (llmax + 3)][iter] * Y2;
					pows[i + 2 * (llmax + 3)][iter] =
						pows[(i - 1) + 2 * (llmax + 3)][iter] * Y3;
				}

				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 4; j++) {
						poly_vgl[i + 5 * j][iter] = 0.;
					}
				}

				poly_vgl[0][iter] = 1.;

				poly_vgl[5][iter] = pows[3][iter];
				poly_vgl[6][iter] = 1.;

				poly_vgl[10][iter] = pows[3 + (llmax + 3)][iter];
				poly_vgl[12][iter] = 1.;

				poly_vgl[15][iter] = pows[3 + 2 * (llmax + 3)][iter];
				poly_vgl[18][iter] = 1.;

				n = 3;
			}

			double xy, yz, xz;
			double da, db, dc, dd;

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

						xy = pows[(a + 2)][iter] *
							 pows[(b + 2) + (llmax + 3)][iter];
						yz = pows[(b + 2) + (llmax + 3)][iter] *
							 pows[(c + 2) + 2 * (llmax + 3)][iter];
						xz = pows[(a + 2)][iter] *
							 pows[(c + 2) + 2 * (llmax + 3)][iter];

						poly_vgl[5 * (n)][iter] =
							xy * pows[c + 2 + 2 * (llmax + 3)][iter];

						xy = dc * xy;
						xz = db * xz;
						yz = da * yz;

						poly_vgl[1 + 5 * n][iter] = pows[a + 1][iter] * yz;
						poly_vgl[2 + 5 * n][iter] =
							pows[b + 1 + (llmax + 3)][iter] * xz;
						poly_vgl[3 + 5 * n][iter] =
							pows[c + 1 + 2 * (llmax + 3)][iter] * xy;

						poly_vgl[4 + 5 * n][iter] =
							(da - 1.) * pows[a][iter] * yz +
							(db - 1.) * pows[b + (llmax + 3)][iter] * xz +
							(dc - 1.) * pows[c + 2 * (llmax + 3)][iter] * xy;

						db -= 1.;
					}
					da -= 1.;
				}
				dd += 1.;
			}
			}
			// End of ao_polynomial computation (now inlined)
			// poly_vgl is now set from here
		
#pragma acc parallel loop collapse(2) async(1) \
	deviceptr(\
	ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,         \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl,    \
	pows, nucl_coord, pows_shared, shell_ang_mom, nucleus_range,               \
	shell_to_nucl)
#pragma omp target teams loop collapse(2) \
	is_device_ptr(\
	ao_value, lstart, ao_index, ao_factor, coord, nucleus_max_ang_mom,         \
	nucleus_index, nucleus_shell_num, shell_vgl, poly_vgl_shared, poly_vgl,    \
	pows, nucl_coord, pows_shared, shell_ang_mom, nucleus_range,               \
	shell_to_nucl)
		for (int iter_new = 0; iter_new < chunk_size / nucl_num; iter_new++) {
			for (int ishell = 0; ishell < shell_num; ishell++) {

				double(*poly_vgl)[chunk_size] =
					(double(*)[chunk_size])poly_vgl_shared;

				int ipoint = iter_new + chunk_size / nucl_num * sub_iter;
				int inucl = shell_to_nucl[ishell];
				int iter = iter_new * nucl_num + inucl;

				if (ipoint >= point_num)
					continue;

				double e_coord_0 = coord[0 * point_num + ipoint];
				double e_coord_1 = coord[1 * point_num + ipoint];
				double e_coord_2 = coord[2 * point_num + ipoint];

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

				int k = ao_index[ishell] - 1;
				int l = shell_ang_mom[ishell];

				for (int il = lstart[l] - 1; il <= lstart[l + 1] - 2; il++) {

					ao_value[k + ipoint * ao_num] =
						poly_vgl[il * 5][iter] *
						shell_vgl[ishell + ipoint * shell_num * 5] *
						ao_factor[k];

					k = k + 1;
				}
			}
		}
#pragma acc wait(1)
	// End of outer compute loop
}
// End of target data region
qmckl_free_device(context, ao_index);
qmckl_free_device(context, poly_vgl_shared);
qmckl_free_device(context, pows_shared);
qmckl_free_device(context, shell_to_nucl);
qmckl_free_device(context, lstart);

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

#pragma acc kernels deviceptr(v, vgl)
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

	// TODO Specify OpenACC clauses manually

#pragma acc kernels copy(shell_max_ptr[:1], prim_max_ptr[:1]) deviceptr(          \
	nucleus_shell_num, nucleus_index, shell_prim_num, prim_num_per_nucleus)
#pragma omp target map(tofrom : shell_max_ptr[:1], prim_max_ptr[:1]) \
	is_device_ptr(nucleus_shell_num, nucleus_index, shell_prim_num,  \
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

	// TODO Manually write OpenACC clauses

#pragma acc kernels deviceptr(expo_expo, expo_index, coef, newcoef,               \
						   nucleus_index, shell_prim_index, nucleus_shell_num, \
						   exponent, coefficient_normalized, shell_prim_num,   \
						   expo_per_nucleus_data, coef_per_nucleus_data,       \
						   prim_num_per_nucleus, newidx)
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

// TODO Manually add OpenACC clauses
#pragma acc kernels deviceptr(nucleus_index, nucleus_prim_index, shell_prim_index)
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

#pragma acc kernels deviceptr(shell_prim_index, shell_prim_num,                   \
						   coefficient_normalized, coefficient, prim_factor,   \
						   shell_factor)
#pragma omp target is_device_ptr(shell_prim_index, shell_prim_num,             \
								 coefficient_normalized, coefficient,          \
								 prim_factor, shell_factor)
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

#pragma acc kernels deviceptr(nucleus_max_ang_mom, nucleus_index,                 \
						   nucleus_shell_num, shell_ang_mom)
#pragma omp target is_device_ptr(nucleus_max_ang_mom, nucleus_index,           \
								 nucleus_shell_num, shell_ang_mom)
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

#pragma acc kernels deviceptr(nucleus_range, nucleus_index, nucleus_shell_num,    \
						   shell_prim_index, shell_prim_num, exponent)
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

	rc = qmckl_finalize_ao_basis_hpc_device(context);

	return rc;
}

qmckl_exit_code_device
qmckl_init_ao_basis_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_init_ao_basis_device", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	ctx->ao_basis.uninitialized = (1 << 14) - 1;

	/* Default values */
	ctx->ao_basis.ao_cartesian = true;

	return QMCKL_SUCCESS_DEVICE;
}

/* Provided check  */

bool qmckl_ao_basis_provided_device(qmckl_context_device context) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->ao_basis.provided;
}

//**********
// PROVIDE
//**********

/* shell_vgl */

qmckl_exit_code_device
qmckl_provide_ao_basis_shell_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_ao_basis_ao_vgl_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_provide_ao_basis_shell_vgl_device",
									 NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.shell_vgl_date) {

		/* Allocate array */
		if (ctx->ao_basis.shell_vgl == NULL) {

			double *shell_vgl = (double *)qmckl_malloc_device(
				context,
				ctx->ao_basis.shell_num * 5 * ctx->point.num * sizeof(double));

			if (shell_vgl == NULL) {
				return qmckl_failwith_device(
					context, QMCKL_ALLOCATION_FAILED_DEVICE,
					"qmckl_ao_basis_shell_vgl_device", NULL);
			}
			ctx->ao_basis.shell_vgl = shell_vgl;
		}

		qmckl_exit_code_device rc;
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
			return qmckl_failwith_device(
				context, QMCKL_FAILURE_DEVICE, "compute_ao_basis_shell_vgl",
				"Not yet implemented for basis type != 'G'");
		}
		if (rc != QMCKL_SUCCESS_DEVICE) {

			return rc;
		}

		ctx->ao_basis.shell_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_vgl */

qmckl_exit_code_device
qmckl_provide_ao_basis_ao_vgl_device(qmckl_context_device context) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_provide_ao_basis_ao_vgl_device",
									 NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	if (!ctx->ao_basis.provided) {
		return qmckl_failwith_device((qmckl_context_device)context,
									 QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_ao_basis_ao_vgl_device", NULL);
	}

	/* Compute if necessary */
	if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

		qmckl_exit_code_device rc;

		/* Allocate array */
		if (ctx->ao_basis.ao_vgl == NULL) {

			double *ao_vgl = (double *)qmckl_malloc_device(
				context,
				ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double));

			if (ao_vgl == NULL) {
				return qmckl_failwith_device(context,
											 QMCKL_ALLOCATION_FAILED_DEVICE,
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
				ctx->point.num, ctx->nucleus.num, ctx->point.coord.data,
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

		if (rc != QMCKL_SUCCESS_DEVICE) {
			return rc;
		}

		ctx->ao_basis.ao_vgl_date = ctx->date;
	}

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// GET
//**********

/* shell_vgl */

qmckl_exit_code_device
qmckl_get_ao_basis_shell_vgl_device(qmckl_context_device context,
									double *const shell_vgl,
									const int64_t size_max) {
	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_vgl_device",
									 NULL);
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ao_basis_shell_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	// int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.shell_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_vgl",
									 "input array too small");
	}
	qmckl_memcpy_D2D(context, shell_vgl, ctx->ao_basis.shell_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_vgl_gaussian */

qmckl_exit_code_device
qmckl_get_ao_basis_ao_vgl_device(qmckl_context_device context,
								 double *const ao_vgl, const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device((qmckl_context_device)context,
									 QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_ao_vgl_device", NULL);
	}

	qmckl_exit_code_device rc;
	rc = qmckl_provide_ao_basis_ao_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_ao_vgl_device", "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_vgl, ctx->ao_basis.ao_vgl,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_vgl_inplace_device(qmckl_context_device context,
										 double *const ao_vgl,
										 const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_ao_vgl_inplace_device", NULL);
	}

	qmckl_exit_code_device rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device((qmckl_context_device)context,
									 QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_ao_vgl_inplace_device",
									 "input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	double *old_array = ctx->ao_basis.ao_vgl;

	ctx->ao_basis.ao_vgl = ao_vgl;

	rc = qmckl_provide_ao_basis_ao_vgl_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->ao_basis.ao_vgl = old_array;

	return QMCKL_SUCCESS_DEVICE;
}

/* ao_value */

qmckl_exit_code_device
qmckl_get_ao_basis_ao_value_device(qmckl_context_device context,
								   double *const ao_value,
								   const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_ao_value_device", NULL);
	}

	qmckl_exit_code_device rc;

	rc = qmckl_provide_ao_basis_ao_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_ao_value_device", "input array too small");
	}

	qmckl_memcpy_D2D(context, ao_value, ctx->ao_basis.ao_value,
					 (size_t)sze * sizeof(double));

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_value_inplace_device(qmckl_context_device context,
										   double *const ao_value,
										   const int64_t size_max) {

	if (qmckl_context_check_device((qmckl_context_device)context) ==
		QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_ao_value_inplace_device", NULL);
	}

	qmckl_exit_code_device rc;

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);
	int device_id = qmckl_get_device_id(context);

	int64_t sze = ctx->ao_basis.ao_num * ctx->point.num;
	if (size_max < sze) {
		return qmckl_failwith_device(
			(qmckl_context_device)context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_ao_value_inplace_device",
			"input array too small");
	}

	rc = qmckl_context_touch_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	double *old_array = ctx->ao_basis.ao_value;

	ctx->ao_basis.ao_value = ao_value;

	rc = qmckl_provide_ao_basis_ao_value_device(context);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	ctx->ao_basis.ao_value = old_array;

	return QMCKL_SUCCESS_DEVICE;
}

//**********
// VARIOUS GETTERS/SETTERS
//**********

qmckl_exit_code_device
qmckl_set_ao_basis_type_device(qmckl_context_device context, char basis_type) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_type_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (basis_type != 'G' && basis_type != 'S') {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_type_device", NULL);
	}

	int32_t mask = 1;

	ctx->ao_basis.type = basis_type;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_shell_num_device(qmckl_context_device context,
									int64_t shell_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (shell_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num <= 0");
	}

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (0L < prim_num && prim_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num > prim_num");
	}

	int32_t mask = 1 << 1;

	ctx->ao_basis.shell_num = shell_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_prim_num_device(qmckl_context_device context,
								   int64_t prim_num) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (prim_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "prim_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "shell_num is not set");
	}

	if (prim_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_prim_num_device",
									 "prim_num < shell_num");
	}

	int32_t mask = 1 << 2;

	ctx->ao_basis.prim_num = prim_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_ao_num_device(qmckl_context_device context, int64_t ao_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_ao_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	if (ao_num <= 0) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "ao_num must be positive");
	}

	int64_t shell_num = ctx->ao_basis.shell_num;
	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "shell_num is not set");
	}

	if (ao_num < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_set_ao_basis_shell_num_device",
									 "ao_num < shell_num");
	}

	int32_t mask = 1 << 12;

	ctx->ao_basis.ao_num = ao_num;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 4;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 "nucl_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 "input array too small");
	}

	if (ctx->ao_basis.nucleus_index != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_index);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_nucleus_index_device", NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_nucleus_index_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_index,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 3;

	int64_t nucl_num = ctx->nucleus.num;

	if (nucl_num <= 0L) {
		return qmckl_failwith_device(
			context, QMCKL_FAILURE_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device",
			"shell_num is not set");
	}

	if (size_max < nucl_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device",
			"input array too small");
	}

	if (ctx->ao_basis.nucleus_shell_num != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.nucleus_shell_num);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_nucleus_shell_num_device",
				NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, nucl_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_ALLOCATION_FAILED_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_num_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, nucleus_shell_num,
					 nucl_num * sizeof(int64_t));

	ctx->ao_basis.nucleus_shell_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 5;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num == 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_ang_mom != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_ang_mom);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_ang_mom_device", NULL);
		}
	}

	int32_t *new_array =
		(int32_t *)qmckl_malloc_device(context, shell_num * sizeof(int32_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_ang_mom,
					 shell_num * sizeof(int32_t));

	ctx->ao_basis.shell_ang_mom = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_nucleus_shell_prim_num_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 6;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_prim_num != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_num);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_prim_num_device", NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_prim_num_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_num,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_num = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_set_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device", NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 7;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(
			context, QMCKL_FAILURE_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device",
			"shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device",
			"input array too small");
	}

	if (ctx->ao_basis.shell_prim_index != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_prim_index);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_prim_index_device",
				NULL);
		}
	}

	int64_t *new_array =
		(int64_t *)qmckl_malloc_device(context, shell_num * sizeof(int64_t));

	if (new_array == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_ALLOCATION_FAILED_DEVICE,
			"qmckl_set_ao_basis_shell_prim_index_device", NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_prim_index,
					 shell_num * sizeof(int64_t));

	ctx->ao_basis.shell_prim_index = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 8;

	int64_t shell_num = ctx->ao_basis.shell_num;

	if (shell_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 "shell_num is not set");
	}

	if (size_max < shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.shell_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.shell_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_shell_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, shell_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_shell_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, shell_factor,
					 shell_num * sizeof(double));

	ctx->ao_basis.shell_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 9;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 "input array too small");
	}

	if (ctx->ao_basis.exponent != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.exponent);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_exponent_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_exponent_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, exponent, prim_num * sizeof(double));

	ctx->ao_basis.exponent = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 10;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 "input array too small");
	}

	if (ctx->ao_basis.coefficient != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.coefficient);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_coefficient_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_coefficient_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, coefficient,
					 prim_num * sizeof(double));

	ctx->ao_basis.coefficient = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	// Accepts an host array and copies it on device

	int32_t mask = 1 << 11;

	int64_t prim_num = ctx->ao_basis.prim_num;

	if (prim_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 "prim_num is not set");
	}

	if (size_max < prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.prim_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.prim_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_prim_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, prim_factor,
					 prim_num * sizeof(double));

	ctx->ao_basis.prim_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_set_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_set_ao_basis_prim_factor_device",
									 NULL);
	}
	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;

	int device_id = qmckl_get_device_id(context);

	int32_t mask = 1 << 13;

	int64_t ao_num = ctx->ao_basis.ao_num;

	if (ao_num <= 0L) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 "ao_num is not set");
	}

	if (size_max < ao_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 "input array too small");
	}

	if (ctx->ao_basis.ao_factor != NULL) {
		qmckl_exit_code_device rc =
			qmckl_free_device(context, ctx->ao_basis.ao_factor);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			return qmckl_failwith_device(
				context, rc, "qmckl_set_ao_basis_ao_factor_device", NULL);
		}
	}

	double *new_array =
		(double *)qmckl_malloc_device(context, ao_num * sizeof(double));

	if (new_array == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_set_ao_basis_ao_factor_device",
									 NULL);
	}

	qmckl_memcpy_D2D(context, new_array, ao_factor, ao_num * sizeof(double));

	ctx->ao_basis.ao_factor = new_array;

	ctx->ao_basis.uninitialized &= ~mask;
	ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
	if (ctx->ao_basis.provided) {
		qmckl_exit_code_device rc_ = qmckl_finalize_ao_basis_device(context);
		if (rc_ != QMCKL_SUCCESS_DEVICE)
			return rc_;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *ao_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_ao_num", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 12;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_ao_num", NULL);
	}

	assert(ctx->ao_basis.ao_num > (int64_t)0);

	*ao_num = ctx->ao_basis.ao_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_shell_num_device(const qmckl_context_device context,
									int64_t *const shell_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_factor", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_num", NULL);
	}

	assert(ctx->ao_basis.shell_num > (int64_t)0);
	*shell_num = ctx->ao_basis.shell_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_prim_num_device(const qmckl_context_device context,
								   int64_t *prim_num) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_prim_num", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 2;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_prim_num", NULL);
	}

	assert(ctx->ao_basis.prim_num > (int64_t)0);

	*prim_num = ctx->ao_basis.prim_num;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_type_device(qmckl_context_device context,
							   char *const basis_type) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_type", NULL);
	}

	qmckl_context_struct_device *const ctx =
		(qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_type", NULL);
	}

	assert(ctx->ao_basis.type != (char)0);

	basis_type[0] = ctx->ao_basis.type;
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_nucleus_shell_num_device(qmckl_context_device context,
											int64_t *nucleus_shell_num,
											int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 3;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(
			context, QMCKL_NOT_PROVIDED_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", NULL);
	}

	if (nucleus_shell_num == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_2_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device", "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_nucleus_shell_num_device",
			"Array too small. Expected nucl_num");
	}

	assert(ctx->ao_basis.nucleus_shell_num != NULL);
	qmckl_memcpy_D2D(context, nucleus_shell_num,
					 ctx->ao_basis.nucleus_shell_num,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_nucleus_index_device(
	qmckl_context_device context, int64_t *nucleus_index, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 4;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 NULL);
	}

	if (nucleus_index == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 "NULL pointer");
	}

	if (size_max < ctx->nucleus.num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_nucleus_index_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.nucleus_index != NULL);
	qmckl_memcpy_D2D(context, nucleus_index, ctx->ao_basis.nucleus_index,
					 (size_t)ctx->nucleus.num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_ang_mom_device(
	qmckl_context_device context, int32_t *shell_ang_mom, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 5;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 NULL);
	}

	if (shell_ang_mom == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_ang_mom_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_ang_mom != NULL);
	qmckl_memcpy_D2D(context, shell_ang_mom, ctx->ao_basis.shell_ang_mom,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int32_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_shell_factor_device(qmckl_context_device context,
									   double *shell_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 8;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 NULL);
	}

	if (shell_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_factor_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_factor != NULL);
	qmckl_memcpy_D2D(context, shell_factor, ctx->ao_basis.shell_factor,
					 (size_t)ctx->ao_basis.shell_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_num_device(
	qmckl_context_device context, int64_t *shell_prim_num, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 6;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 NULL);
	}

	if (shell_prim_num == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_shell_prim_num_device",
									 "Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_num != NULL);
	qmckl_memcpy_D2D(context, shell_prim_num, ctx->ao_basis.shell_prim_num,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_get_ao_basis_shell_prim_index_device(
	qmckl_context_device context, int64_t *shell_prim_index, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_CONTEXT_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 7;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(
			context, QMCKL_NOT_PROVIDED_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", NULL);
	}

	if (shell_prim_index == NULL) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_2_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device", "NULL pointer");
	}

	if (size_max < ctx->ao_basis.shell_num) {
		return qmckl_failwith_device(
			context, QMCKL_INVALID_ARG_3_DEVICE,
			"qmckl_get_ao_basis_shell_prim_index_device",
			"Array too small. Expected shell_num");
	}

	assert(ctx->ao_basis.shell_prim_index != NULL);
	qmckl_memcpy_D2D(context, shell_prim_index, ctx->ao_basis.shell_prim_index,
					 (size_t)ctx->ao_basis.shell_num * sizeof(int64_t));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_exponent_device(qmckl_context_device context,
								   double *exponent, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 9;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 NULL);
	}

	if (exponent == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_exponent_device",
									 "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.exponent != NULL);
	qmckl_memcpy_D2D(context, exponent, ctx->ao_basis.exponent,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_coefficient_device(qmckl_context_device context,
									  double *coefficient, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 10;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 NULL);
	}

	if (coefficient == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_coefficient_device",
									 "Array too small. Expected prim_num");
	}
	assert(ctx->ao_basis.coefficient != NULL);
	qmckl_memcpy_D2D(context, coefficient, ctx->ao_basis.coefficient,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_prim_factor_device(qmckl_context_device context,
									  double *prim_factor, int64_t size_max) {

	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 11;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 NULL);
	}

	if (prim_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.prim_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_prim_factor_device",
									 "Array too small. Expected prim_num");
	}

	assert(ctx->ao_basis.prim_factor != NULL);
	qmckl_memcpy_D2D(context, prim_factor, ctx->ao_basis.prim_factor,
					 (size_t)ctx->ao_basis.prim_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_get_ao_basis_ao_factor_device(qmckl_context_device context,
									double *ao_factor, int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return qmckl_failwith_device(context, QMCKL_INVALID_CONTEXT_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 NULL);
	}

	qmckl_context_struct_device *ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 13;

	if ((ctx->ao_basis.uninitialized & mask) != 0) {
		return qmckl_failwith_device(context, QMCKL_NOT_PROVIDED_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 NULL);
	}

	if (ao_factor == NULL) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 "NULL pointer");
	}

	if (size_max < ctx->ao_basis.ao_num) {
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_3_DEVICE,
									 "qmckl_get_ao_basis_ao_factor_device",
									 "Array too small. Expected ao_num");
	}

	assert(ctx->ao_basis.ao_factor != NULL);
	qmckl_memcpy_D2D(context, ao_factor, ctx->ao_basis.ao_factor,
					 ctx->ao_basis.ao_num * sizeof(double));
	return QMCKL_SUCCESS_DEVICE;
}
