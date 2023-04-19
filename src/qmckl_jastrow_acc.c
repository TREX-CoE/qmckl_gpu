#include "../include/qmckl_jastrow.h"

//**********
// COMPUTES
//**********

// Finalize computes
qmckl_exit_code_device qmckl_compute_jastrow_champ_asymp_jasa(
	const qmckl_context_device context, const int64_t aord_num,
	const int64_t type_nucl_num, const double *a_vector,
	const double *rescale_factor_en, double *const asymp_jasa) {

	int i, j, p;
	float kappa_inv, x, asym_one;
	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (aord_num < 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	for (int i = 0; i < type_nucl_num; i++) {

		kappa_inv = 1.0 / rescale_factor_en[i];

		asymp_jasa[i] = a_vector[0 + i * (aord_num + 1)] * kappa_inv /
						(1.0 + a_vector[1 + i * (aord_num + 1)] * kappa_inv);

		x = kappa_inv;
		for (int p = 1; p < aord_num; p++) {
			x = x * kappa_inv;
			asymp_jasa[i] =
				asymp_jasa[i] + a_vector[p + 1 + i * (aord_num + 1)] * x;
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_jastrow_asymp_jasb_device(
	const qmckl_context_device context, const int64_t bord_num,
	const double *b_vector, const double rescale_factor_ee,
	double *const asymp_jasb) {

	double asym_one, x;
	double kappa_inv = 1.0 / rescale_factor_ee;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (bord_num < 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	asym_one = b_vector[0] * kappa_inv / (1.0 + b_vector[1] * kappa_inv);

	asymp_jasb[0] = asym_one;
	asymp_jasb[1] = 0.5 * asym_one;

	for (int i = 0; i < 2; i++) {
		x = kappa_inv;
		for (int p = 1; p < bord_num; p++) {
			x = x * kappa_inv;
			asymp_jasb[i] = asymp_jasb[i] + b_vector[p + 1] * x;
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Total Jastrow
qmckl_exit_code_device
qmckl_compute_jastrow_value_device(const qmckl_context_device context,
								   const int64_t walk_num, const double *f_ee,
								   const double *f_en, const double *f_een,
								   double *const value) {
	// TODO
	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron component
qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *asymp_jasb, double *const factor_ee) {
	// TODO
	return QMCKL_SUCCESS_DEVICE;
}

// Electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_en_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *asymp_jasa, double *const factor_en) {
	// TODO
	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_een_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *een_rescaled_e,
	const double *een_rescaled_n, double *const factor_een) {
	// TODO
	return QMCKL_SUCCESS_DEVICE;
}

// Distances
qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled) {

  int k;

  qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

  if (context == QMCKL_NULL_CONTEXT_DEVICE) {
     info = QMCKL_INVALID_CONTEXT_DEVICE;
     return info;
  }

  if (elec_num <= 0) {
     info = QMCKL_INVALID_ARG_2_DEVICE;
     return info;
  }

  if (walk_num <= 0) {
     info = QMCKL_INVALID_ARG_3_DEVICE;
     return info;
  }

  for(int k=0; k<walk_num; k++) {
     // TODO
     info = qmckl_distance_rescaled_device(context, 'T', 'T', elec_num, elec_num,
          coord[k * elec_num], elec_num * walk_num,
          coord[k * elec_num], elec_num * walk_num,
          ee_distance_rescaled[k * elec_num * elec_num], elec_num, rescale_factor_ee);
     if (info != QMCKL_SUCCESS_DEVICE) {
        return info;
     }
  }

}

qmckl_exit_code_device qmckl_compute_en_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled);

qmckl_exit_code_device qmckl_compute_een_rescaled_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t cord_num,
	const double rescale_factor_ee, const double *ee_distance,
	double *const een_rescaled_e) {
	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (walk_num <= 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3_DEVICE;
	}

	if (cord_num < 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

	// Prepare table of exponentiated distances raised to appropriate power
	// init

	memset(een_rescaled_e, 0,
		   walk_num * (cord_num + 1) * elec_num * elec_num * sizeof(double));

	const size_t elec_pairs = (size_t)(elec_num * (elec_num - 1)) / 2;
	const size_t len_een_ij = (size_t)elec_pairs * (cord_num + 1);

	// number of elements for the een_rescaled_e_ij[N_e*(N_e-1)/2][cord+1]
	// probably in C is better [cord+1, Ne*(Ne-1)/2]
	// elec_pairs = (elec_num * (elec_num - 1)) / 2;
	// len_een_ij = elec_pairs * (cord_num + 1);
	const size_t e2 = elec_num * elec_num;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
	for (size_t nw = 0; nw < (size_t)walk_num; ++nw) {

		double een_rescaled_e_ij[len_een_ij];

		memset(&(een_rescaled_e_ij[0]), 0, len_een_ij * sizeof(double));
		for (size_t kk = 0; kk < elec_pairs; ++kk) {
			een_rescaled_e_ij[kk] = 1.0;
		}

		size_t kk = 0;
		for (size_t i = 0; i < (size_t)elec_num; ++i) {
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
			for (size_t j = 0; j < i; ++j) {
				een_rescaled_e_ij[j + kk + elec_pairs] =
					-rescale_factor_ee *
					ee_distance[j + i * elec_num + nw * e2];
			}
			kk += i;
		}

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
		for (size_t k = elec_pairs; k < 2 * elec_pairs; ++k) {
			een_rescaled_e_ij[k] = exp(een_rescaled_e_ij[k]);
		}

		for (size_t l = 2; l < (size_t)(cord_num + 1); ++l) {
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
			for (size_t k = 0; k < elec_pairs; ++k) {
				// een_rescaled_e_ij(k, l + 1) = een_rescaled_e_ij(k, l + 1 - 1)
				// * een_rescaled_e_ij(k, 2)
				een_rescaled_e_ij[k + l * elec_pairs] =
					een_rescaled_e_ij[k + (l - 1) * elec_pairs] *
					een_rescaled_e_ij[k + elec_pairs];
			}
		}

		double *const een_rescaled_e_ =
			&(een_rescaled_e[nw * (cord_num + 1) * e2]);
		// prepare the actual een table
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
		for (size_t i = 0; i < e2; ++i) {
			een_rescaled_e_[i] = 1.0;
		}

		for (size_t l = 1; l < (size_t)(cord_num + 1); ++l) {
			double *x = een_rescaled_e_ij + l * elec_pairs;
			double *const een_rescaled_e__ = &(een_rescaled_e_[l * e2]);
			double *een_rescaled_e_i = een_rescaled_e__;
			for (size_t i = 0; i < (size_t)elec_num; ++i) {
				for (size_t j = 0; j < i; ++j) {
					een_rescaled_e_i[j] = *x;
					een_rescaled_e__[i + j * elec_num] = *x;
					x += 1;
				}
				een_rescaled_e_i += elec_num;
			}
		}

		double *const x0 = &(een_rescaled_e[nw * e2 * (cord_num + 1)]);
		for (size_t l = 0; l < (size_t)(cord_num + 1); ++l) {
			double *x1 = &(x0[l * e2]);
			for (size_t j = 0; j < (size_t)elec_num; ++j) {
				*x1 = 0.0;
				x1 += 1 + elec_num;
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_een_rescaled_n_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const int64_t cord_num,
	const double *rescale_factor_en, const double *en_distance,
	double *const een_rescaled_n) {
	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (walk_num <= 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3_DEVICE;
	}

	if (nucl_num <= 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

	if (cord_num < 0) {
		return QMCKL_INVALID_ARG_5_DEVICE;
	}

	// Prepare table of exponentiated distances raised to appropriate power
	for (int i = 0; i < (walk_num * (cord_num + 1) * nucl_num * elec_num);
		 ++i) {
		een_rescaled_n[i] = 1.0;
	}

	for (int nw = 0; nw < walk_num; ++nw) {
		for (int a = 0; a < nucl_num; ++a) {
			for (int i = 0; i < elec_num; ++i) {
				een_rescaled_n[i + a * elec_num +
							   nw * elec_num * nucl_num * (cord_num + 1)] = 1.0;
				een_rescaled_n[i + a * elec_num + elec_num * nucl_num +
							   nw * elec_num * nucl_num * (cord_num + 1)] =
					exp(-rescale_factor_en[type_nucl_vector[a]] *
						en_distance[a + i * nucl_num +
									nw * elec_num * nucl_num]);
			}
		}

		for (int l = 2; l < (cord_num + 1); ++l) {
			for (int a = 0; a < nucl_num; ++a) {
				for (int i = 0; i < elec_num; ++i) {
					een_rescaled_n[i + a * elec_num + l * elec_num * nucl_num +
								   nw * elec_num * nucl_num * (cord_num + 1)] =
						een_rescaled_n[i + a * elec_num +
									   (l - 1) * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)] *
						een_rescaled_n[i + a * elec_num + elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)];
				}
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_c_vector_full_device(
	const qmckl_context_device context, const int64_t nucl_num,
	const int64_t dim_c_vector, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const double *c_vector,
	double *const c_vector_full) {
	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (nucl_num <= 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (type_nucl_num <= 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

	if (dim_c_vector < 0) {
		return QMCKL_INVALID_ARG_5_DEVICE;
	}

	for (int i = 0; i < dim_c_vector; ++i) {
		for (int a = 0; a < nucl_num; ++a) {
			c_vector_full[a + i * nucl_num] =
				c_vector[(type_nucl_vector[a] - 1) + i * type_nucl_num];
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_lkpm_combined_index_device(
	const qmckl_context_device context, const int64_t cord_num,
	const int64_t dim_c_vector, int64_t *const lkpm_combined_index) {

	double x;
	int i, a, k, l, kk, p, lmax, m;

	qmckl_context_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (cord_num < 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (dim_c_vector < 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	kk = 0;
	for (p = 1; p < cord_num; p++) {
		for (k = p - 1 - 2; k >= 0; k--) {
			if (k != 0) {
				lmax = p - k;
			} else {
				lmax = p - k - 2;
			}
			for (l = lmax - 1; lmax >= 0; l--) {
				if ((p - k - l) & 1 == 1)
					continue;
				m = (p - k - l) / 2;
				kk = kk + 1;
				lkpm_combined_index[kk] = l;
				lkpm_combined_index[kk + dim_c_vector] = k;
				lkpm_combined_index[kk + 2 * dim_c_vector] = p;
				lkpm_combined_index[kk + 3 * dim_c_vector] = m;
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_compute_tmp_c_device(const qmckl_context_device context,
						   const int64_t cord_num, const int64_t elec_num,
						   const int64_t nucl_num, const int64_t walk_num,
						   const double *een_rescaled_e,
						   const double *een_rescaled_n, double *const tmp_c) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (cord_num < 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3_DEVICE;
	}

	if (nucl_num <= 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

	if (walk_num <= 0) {
		return QMCKL_INVALID_ARG_5_DEVICE;
	}

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	const char TransA = 'N';
	const char TransB = 'N';
	const double alpha = 1.0;
	const double beta = 0.0;

	const int64_t M = elec_num;
	const int64_t N = nucl_num * (cord_num + 1);
	const int64_t K = elec_num;

	const int64_t LDA = elec_num;
	const int64_t LDB = elec_num;
	const int64_t LDC = elec_num;

	const int64_t af = elec_num * elec_num;
	const int64_t bf = elec_num * nucl_num * (cord_num + 1);
	const int64_t cf = bf;

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
	for (int64_t nw = 0; nw < walk_num; ++nw) {
		for (int64_t i = 0; i < cord_num; ++i) {
			/* TODO Replace by a BLAS call or write a temporary manual DGEMM
			info =
				qmckl_dgemm(context, TransA, TransB, M, N, K, alpha,
							&(een_rescaled_e[af * (i + nw * (cord_num + 1))]),
							LDA, &(een_rescaled_n[bf * nw]), LDB, beta,
							&(tmp_c[cf * (i + nw * cord_num)]), LDC);
			 */
		}
	}

	return info;
}
