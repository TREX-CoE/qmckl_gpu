#include "../include/qmckl_jastrow.h"

//**********
// COMPUTES
//**********

// Finalize computes
qmckl_exit_code_device qmckl_compute_jastrow_asymp_jasa_device(
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

#pragma acc kernels deviceptr(asymp_jasa, a_vector, rescale_factor_en)
	{
		for (int i = 0; i < type_nucl_num; i++) {

			kappa_inv = 1.0 / rescale_factor_en[i];

			asymp_jasa[i] =
				a_vector[0 + i * (aord_num + 1)] * kappa_inv /
				(1.0 + a_vector[1 + i * (aord_num + 1)] * kappa_inv);

			x = kappa_inv;
			for (int p = 1; p < aord_num; p++) {
				x = x * kappa_inv;
				asymp_jasa[i] =
					asymp_jasa[i] + a_vector[p + 1 + i * (aord_num + 1)] * x;
			}
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

#pragma acc kernels deviceptr(asymp_jasb, b_vector)
	{
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
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Total Jastrow
qmckl_exit_code_device
qmckl_compute_jastrow_value_device(const qmckl_context_device context,
								   const int64_t walk_num, const double *f_ee,
								   const double *f_en, const double *f_een,
								   double *const value) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_INVALID_CONTEXT_DEVICE;
	if (walk_num <= 0)
		return QMCKL_INVALID_ARG_2_DEVICE;
	if (f_ee == NULL)
		return QMCKL_INVALID_ARG_3_DEVICE;
	if (f_en == NULL)
		return QMCKL_INVALID_ARG_4_DEVICE;
	if (f_een == NULL)
		return QMCKL_INVALID_ARG_5_DEVICE;
	if (value == NULL)
		return QMCKL_INVALID_ARG_6_DEVICE;

#pragma acc kernels deviceptr(value, f_ee, f_en, f_een)
	{
		for (int64_t i = 0; i < walk_num; ++i) {
			value[i] = exp(f_ee[i] + f_en[i] + f_een[i]);
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_jastrow_gl_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const double *value, const double *gl_ee,
	const double *gl_en, const double *gl_een, double *const gl) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_INVALID_CONTEXT_DEVICE;
	if (walk_num <= 0)
		return QMCKL_INVALID_ARG_2_DEVICE;
	if (elec_num <= 0)
		return QMCKL_INVALID_ARG_3_DEVICE;
	if (value == NULL)
		return QMCKL_INVALID_ARG_4_DEVICE;
	if (gl_ee == NULL)
		return QMCKL_INVALID_ARG_5_DEVICE;
	if (gl_en == NULL)
		return QMCKL_INVALID_ARG_6_DEVICE;
	if (gl_een == NULL)
		return QMCKL_INVALID_ARG_7_DEVICE;
	if (gl == NULL)
		return QMCKL_INVALID_ARG_8_DEVICE;

#pragma acc kernels deviceptr(value, gl_ee, gl_en, gl_een, gl)
	{
		for (int64_t k = 0; k < walk_num; ++k) {
			for (int64_t j = 0; j < 4; ++j) {
				for (int64_t i = 0; i < elec_num; ++i) {
					gl[i + elec_num * (j + k * 4)] =
						gl_ee[i + elec_num * (j + k * 4)] +
						gl_en[i + elec_num * (j + k * 4)] +
						gl_een[i + elec_num * (j + k * 4)];
				}
			}
			for (int64_t i = 0; i < elec_num; ++i) {
				gl[i + elec_num * (3 + walk_num * 4)] +=
					gl_ee[i + elec_num * (0 + k * 4)] *
						gl_ee[i + elec_num * (0 + k * 4)] +
					gl_ee[i + elec_num * (1 + k * 4)] *
						gl_ee[i + elec_num * (1 + k * 4)] +
					gl_ee[i + elec_num * (2 + k * 4)] *
						gl_ee[i + elec_num * (2 + k * 4)];
			}
			for (int64_t j = 0; j < 4; ++j) {
				for (int64_t i = 0; i < elec_num; ++i) {
					gl[i + elec_num * (j + k * 4)] *= value[k];
				}
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/electron component
qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *asymp_jasb, double *const factor_ee) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (walk_num <= 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3_DEVICE;
	}

	if (bord_num < 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

#pragma acc kernels deviceptr(ee_distance_rescaled, factor_ee, b_vector)
	{
		for (int nw = 0; nw < walk_num; ++nw) {
			factor_ee[nw] = 0.0; // put init array here.
			size_t ishift = nw * elec_num * elec_num;
			for (int i = 0; i < elec_num; ++i) {
				for (int j = 0; j < i; ++j) {
					double x = ee_distance_rescaled[j + i * elec_num + ishift];
					const double x1 = x;
					double power_ser = 0.0;
					double spin_fact = 1.0;
					int ipar = 0; // index of asymp_jasb

					for (int p = 1; p < bord_num; ++p) {
						x = x * x1;
						power_ser += b_vector[p + 1] * x;
					}

					if (i < up_num || j >= up_num) {
						spin_fact = 0.5;
						ipar = 1;
					}

					factor_ee[nw] += spin_fact * b_vector[0] * x1 /
										 (1.0 + b_vector[1] * x1) -
									 asymp_jasb[ipar] + power_ser;
				}
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device qmckl_compute_jastrow_factor_ee_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t up_num, const int64_t bord_num,
	const double *b_vector, const double *ee_distance_rescaled,
	const double *ee_distance_rescaled_deriv_e,
	double *const factor_ee_deriv_e) {

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (walk_num <= 0) {
		return QMCKL_INVALID_ARG_2_DEVICE;
	}

	if (elec_num <= 0) {
		return QMCKL_INVALID_ARG_3_DEVICE;
	}

	if (bord_num < 0) {
		return QMCKL_INVALID_ARG_4_DEVICE;
	}

#pragma acc kernels deviceptr(b_vector, ee_distance_rescaled,                  \
								  ee_distance_rescaled_deriv_e,                \
								  factor_ee_deriv_e)
	{

		for (int nw = 0; nw < walk_num; ++nw) {
			for (int ii = 0; ii < 4; ++ii) {
				for (int j = 0; j < elec_num; ++j) {
					factor_ee_deriv_e[j + ii * elec_num + nw * elec_num * 4] =
						0.0;
				}
			}
		}

		const double third = 1.0 / 3.0;

		for (int nw = 0; nw < walk_num; ++nw) {
			for (int i = 0; i < elec_num; ++i) {
				for (int j = 0; j < elec_num; ++j) {
					const double x0 =
						ee_distance_rescaled[j + i * elec_num +
											 nw * elec_num * elec_num];
					if (fabs(x0) < 1.0e-18)
						continue;
					double spin_fact = 1.0;
					const double den = 1.0 + b_vector[1] * x0;
					const double invden = 1.0 / den;
					const double invden2 = invden * invden;
					const double invden3 = invden2 * invden;
					const double xinv = 1.0 / (x0 + 1.0e-18);

					double dx[4];
					dx[0] = ee_distance_rescaled_deriv_e[0 + j * 4 +
														 i * 4 * elec_num +
														 nw * 4 * elec_num *
															 elec_num];
					dx[1] = ee_distance_rescaled_deriv_e[1 + j * 4 +
														 i * 4 * elec_num +
														 nw * 4 * elec_num *
															 elec_num];
					dx[2] = ee_distance_rescaled_deriv_e[2 + j * 4 +
														 i * 4 * elec_num +
														 nw * 4 * elec_num *
															 elec_num];
					dx[3] = ee_distance_rescaled_deriv_e[3 + j * 4 +
														 i * 4 * elec_num +
														 nw * 4 * elec_num *
															 elec_num];

					if ((i <= (up_num - 1) && j <= (up_num - 1)) ||
						(i > (up_num - 1) && j > (up_num - 1))) {
						spin_fact = 0.5;
					}

					double lap1 = 0.0;
					double lap2 = 0.0;
					double lap3 = 0.0;
					double pow_ser_g[3] = {0., 0., 0.};
					for (int ii = 0; ii < 3; ++ii) {
						double x = x0;
						if (fabs(x) < 1.0e-18)
							continue;
						for (int p = 2; p < bord_num + 1; ++p) {
							const double y = p * b_vector[(p - 1) + 1] * x;
							pow_ser_g[ii] = pow_ser_g[ii] + y * dx[ii];
							lap1 = lap1 + (p - 1) * y * xinv * dx[ii] * dx[ii];
							lap2 = lap2 + y;
							x = x *
								ee_distance_rescaled[j + i * elec_num +
													 nw * elec_num * elec_num];
						}

						lap3 = lap3 - 2.0 * b_vector[1] * dx[ii] * dx[ii];

						factor_ee_deriv_e[i + ii * elec_num +
										  nw * elec_num * 4] +=
							+spin_fact * b_vector[0] * dx[ii] * invden2 +
							pow_ser_g[ii];
					}

					int ii = 3;
					lap2 = lap2 * dx[ii] * third;
					lap3 = lap3 + den * dx[ii];
					lap3 = lap3 * (spin_fact * b_vector[0] * invden3);
					factor_ee_deriv_e[i + ii * elec_num + nw * elec_num * 4] +=
						lap1 + lap2 + lap3;
				}
			}
		}
	}

	return QMCKL_SUCCESS_DEVICE;
}

// Electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_en_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *asymp_jasa, double *const factor_en) {

	int i, a, p, nw;
	double x, power_ser;
	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (type_nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (aord_num < 0) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

#pragma acc kernels deviceptr(type_nucl_vector, a_vector,                      \
								  en_distance_rescaled, asymp_jasa)
	{
		for (nw = 0; nw < walk_num; nw++) {
			factor_en[nw] = 0.0;
			for (a = 0; a < nucl_num; a++) {
				for (i = 0; i < elec_num; i++) {
					x = en_distance_rescaled[i + a * elec_num +
											 nw * elec_num * nucl_num];

					factor_en[nw] =
						factor_en[nw] +
						a_vector[0 + type_nucl_vector[a] * (aord_num + 1)] * x /
							(1.0 + a_vector[1 + type_nucl_vector[a] *
													(aord_num + 1)] *
									   x) -
						asymp_jasa[type_nucl_vector[a]];

					for (p = 2; p < aord_num; p++) {
						x = x * en_distance_rescaled[i + a * elec_num +
													 nw * elec_num * nucl_num];
						factor_en[nw] =
							factor_en[nw] +
							a_vector[p + type_nucl_vector[a] * (aord_num + 1)] *
								x;
					}
				}
			}
		}
	}

	return info;
}

qmckl_exit_code_device qmckl_compute_jastrow_factor_en_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *en_distance_rescaled_deriv_e,
	double *const factor_en_deriv_e) {

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (aord_num < 0) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

	int i, a, p, ipar, nw, ii;
	double x, den, invden, invden2, invden3, xinv;
	double y, lap1, lap2, lap3, third;

	double power_ser_g[3];
	double dx[4];

	for (int i = 0; i < elec_num * 4 * walk_num; i++)
		factor_en_deriv_e[i] = 0.0;
	third = 1.0 / 3.0;

	for (nw = 0; nw < walk_num; nw++) {
		for (a = 0; a < nucl_num; a++) {
			for (i = 0; i < elec_num; i++) {
				x = en_distance_rescaled[i + a * elec_num +
										 nw * elec_num * nucl_num];
				if (abs(x) < 1.0e-18)
					continue;
				power_ser_g[0] = 0.0;
				power_ser_g[1] = 0.0;
				power_ser_g[2] = 0.0;
				den = 1.0 +
					  a_vector[1 + type_nucl_vector[a] * (aord_num + 1)] * x;
				invden = 1.0 / den;
				invden2 = invden * invden;
				invden3 = invden2 * invden;
				xinv = 1.0 / x;

				for (ii = 0; ii < 4; ii++) {
					dx[ii] = en_distance_rescaled_deriv_e[ii + i * 4 +
														  a * 4 * elec_num +
														  nw * 4 * elec_num *
															  nucl_num];
				}

				lap1 = 0.0;
				lap2 = 0.0;
				lap3 = 0.0;
				for (ii = 0; ii < 3; ii++) {
					x = en_distance_rescaled[i + a * elec_num +
											 nw * elec_num * nucl_num];
					for (p = 1; p < aord_num; p++) {
						y = p *
							a_vector[(p + 1) +
									 type_nucl_vector[a] * (aord_num + 1)] *
							x;
						power_ser_g[ii] = power_ser_g[ii] + y * dx[ii];
						lap1 = lap1 + p * y * xinv * dx[ii] * dx[ii];
						lap2 = lap2 + y;
						x = x * en_distance_rescaled[i + a * elec_num +
													 nw * elec_num * nucl_num];
					}

					lap3 = lap3 - 2.0 * a_vector[2 + type_nucl_vector[a]] *
									  dx[ii] * dx[ii];

					factor_en_deriv_e[i + ii * elec_num + nw * elec_num * 4] =
						factor_en_deriv_e[i + ii * elec_num +
										  nw * elec_num * 4] +
						a_vector[0 + type_nucl_vector[a] * (aord_num + 1)] *
							dx[ii] * invden2 +
						power_ser_g[ii];
				}

				ii = 4;
				lap2 = lap2 * dx[ii] * third;
				lap3 = lap3 + den * dx[ii];
				lap3 = lap3 * a_vector[0 + type_nucl_vector[a]] * invden3;
				factor_en_deriv_e[i + ii * elec_num + nw * elec_num * 4] =
					factor_en_deriv_e[i + ii * elec_num + nw * elec_num * 4] +
					lap1 + lap2 + lap3;
			}
		}
	}
}

qmckl_exit_code_device qmckl_compute_en_distance_rescaled_deriv_e_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled_deriv_e) {
	int i, k;
	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

	double coord[3];
	for (int i = 0; i < nucl_num; i++) {
		coord[0] = nucl_coord[i + 0];
		coord[1] = nucl_coord[i + 1];
		coord[2] = nucl_coord[i + 2];
		for (k = 0; k < walk_num; k++) {
			info = qmckl_distance_rescaled_deriv_e_device(
				context, 'T', 'T', elec_num, 1, elec_coord + (k * elec_num),
				elec_num * walk_num, coord, 1,
				en_distance_rescaled_deriv_e + (0 + 0 * 4 + i * 4 * elec_num +
												k * 4 * elec_num * nucl_num),
				elec_num, rescale_factor_en[type_nucl_vector[i]]);
			if (info != QMCKL_SUCCESS_DEVICE) {
				return info;
			}
		}
	}
}

qmckl_exit_code_device qmckl_compute_jastrow_champ_factor_en_deriv_e(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	const int64_t *type_nucl_vector, const int64_t aord_num,
	const double *a_vector, const double *en_distance_rescaled,
	const double *en_distance_rescaled_deriv_e,
	double *const factor_en_deriv_e) {

	int i, a, p, ipar, nw, ii;
	double x, den, invden, invden2, invden3, xinv;
	double y, lap1, lap2, lap3, third;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (aord_num < 0) {
		info = QMCKL_INVALID_ARG_7_DEVICE;
		return info;
	}

#pragma acc kernels deviceptr(                                                 \
		type_nucl_vector, a_vector, en_distance_rescaled,                      \
			en_distance_rescaled_deriv_e, factor_en_deriv_e)
	{
		for (i = 0; i < elec_num * 4 * walk_num; i++)
			factor_en_deriv_e[i] = 0.0;
		third = 1.0 / 3.0;

		double power_ser_g[3];
		double dx[4];

		for (nw = 0; nw < walk_num; nw++) {
			for (a = 0; a < nucl_num; a++) {
				for (i = 0; i < elec_num; i++) {
					x = en_distance_rescaled[i + a * elec_num +
											 nw * elec_num * nucl_num];
					if (abs(x) < 1.0e-18)
						continue;
					power_ser_g[0] = 0.0;
					power_ser_g[1] = 0.0;
					power_ser_g[2] = 0.0;
					den =
						1.0 +
						a_vector[1 + type_nucl_vector[a] * (aord_num + 1)] * x;
					invden = 1.0 / den;
					invden2 = invden * invden;
					invden3 = invden2 * invden;
					xinv = 1.0 / x;

					for (ii = 0; ii < 4; ii++) {
						dx[ii] =
							en_distance_rescaled_deriv_e[ii + i * 4 +
														 a * 4 * elec_num +
														 nw * 4 * elec_num *
															 nucl_num];
					}

					lap1 = 0.0;
					lap2 = 0.0;
					lap3 = 0.0;
					for (ii = 0; ii < 3; ii++) {
						x = en_distance_rescaled[i + a * elec_num +
												 nw * elec_num * nucl_num];
						for (int p = 1; p < aord_num; p++) {
							y = p *
								a_vector[p +
										 type_nucl_vector[a] * (aord_num + 1)] *
								x;
							power_ser_g[ii] = power_ser_g[ii] + y * dx[ii];
							lap1 = lap1 + (p - 1) * y * xinv * dx[ii] * dx[ii];
							lap2 = lap2 + y;
							x = x *
								en_distance_rescaled[i + a * elec_num +
													 nw * elec_num * nucl_num];
						}

						lap3 = lap3 - 2.0 *
										  a_vector[1 + type_nucl_vector[a] *
														   (aord_num + 1)] *
										  dx[ii] * dx[ii];

						factor_en_deriv_e[i + ii * elec_num +
										  nw * elec_num * 4] =
							factor_en_deriv_e[i + ii * elec_num * +nw *
													  elec_num * 4] +
							a_vector[0 + type_nucl_vector[a] * (aord_num + 1)] *
								dx[ii] * invden2 +
							power_ser_g[ii];
					}

					ii = 3;
					lap2 = lap2 * dx[ii] * third;
					lap3 = lap3 + den * dx[ii];
					lap3 = lap3 *
						   a_vector[0 + type_nucl_vector[a] * (aord_num + 1)] *
						   invden3;
					factor_en_deriv_e[i + ii * elec_num + nw * elec_num * 4] =
						factor_en_deriv_e[i + ii * elec_num +
										  nw * elec_num * 4] +
						lap1 + lap2 + lap3;
				}
			}
		}
	}

	return info;
}

// Electron/electron/nucleus component
qmckl_exit_code_device qmckl_compute_jastrow_factor_een_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *tmp_c,
	const double *een_rescaled_n, double *const factor_een) {

	int i, a, j, l, k, p, m, n, nw;
	double accu, accu2, cn;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (cord_num < 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

#pragma acc kernels deviceptr(c_vector_full, lkpm_combined_index, tmp_c,       \
								  een_rescaled_n)
	{
		for (nw = 0; nw < walk_num; nw++) {
			factor_een[nw] = 0.0;
			for (n = 0; n < dim_c_vector; n++) {
				l = lkpm_combined_index[n];
				k = lkpm_combined_index[n + dim_c_vector];
				p = lkpm_combined_index[n + 2 * dim_c_vector];
				m = lkpm_combined_index[n + 3 * dim_c_vector];

				for (a = 0; a < nucl_num; a++) {
					cn = c_vector_full[a + n * nucl_num];
					if (cn == 0.0)
						continue;

					accu = 0.0;
					for (int j = 0; j < elec_num; j++) {
						accu =
							accu +
							een_rescaled_n[j + a * elec_num +
										   m * elec_num * nucl_num +
										   nw * elec_num * nucl_num *
											   (cord_num + 1)] *
								tmp_c[j + a * elec_num +
									  (m + l) * elec_num * nucl_num +
									  k * elec_num * nucl_num * (cord_num + 1) +
									  nw * elec_num * nucl_num *
										  (cord_num + 1) * cord_num];
					}
					factor_een[nw] = factor_een[nw] + accu * cn;
				}
			}
		}
	}

	return info;
}

qmckl_exit_code_device qmckl_compute_jastrow_factor_een_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t cord_num,
	const int64_t dim_c_vector, const double *c_vector_full,
	const int64_t *lkpm_combined_index, const double *tmp_c,
	const double *dtmp_c, const double *een_rescaled_n,
	const double *een_rescaled_n_deriv_e, double *const factor_een_deriv_e) {

	int64_t info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_INVALID_CONTEXT_DEVICE;
	if (walk_num <= 0)
		return QMCKL_INVALID_ARG_2_DEVICE;
	if (elec_num <= 0)
		return QMCKL_INVALID_ARG_3_DEVICE;
	if (nucl_num <= 0)
		return QMCKL_INVALID_ARG_4_DEVICE;
	if (cord_num < 0)
		return QMCKL_INVALID_ARG_5_DEVICE;

	double *tmp3 = qmckl_malloc_device(context, elec_num * sizeof(double));
#pragma acc kernels deviceptr(tmp3, c_vector_full, lkpm_combined_index, tmp_c, \
								  dtmp_c, een_rescaled_n,                      \
								  een_rescaled_n_deriv_e, factor_een_deriv_e)
	{

		for (int i = 0; i < elec_num * 4 * walk_num; i++) {
			factor_een_deriv_e[i] = 0.;
		}

		const size_t elec_num2 = elec_num << 1;
		const size_t elec_num3 = elec_num * 3;

		for (size_t nw = 0; nw < (size_t)walk_num; ++nw) {
			double *const restrict factor_een_deriv_e_0nw =
				&(factor_een_deriv_e[elec_num * 4 * nw]);
			for (size_t n = 0; n < (size_t)dim_c_vector; ++n) {
				const size_t l = lkpm_combined_index[n];
				const size_t k = lkpm_combined_index[n + dim_c_vector];
				const size_t m = lkpm_combined_index[n + 3 * dim_c_vector];

				const size_t en = elec_num * nucl_num;
				const size_t len = l * en;
				const size_t len4 = len << 2;
				const size_t cn = cord_num * nw;
				const size_t c1 = cord_num + 1;
				const size_t addr0 = en * (m + c1 * (k + cn));
				const size_t addr1 = en * (m + cn);

				const double *restrict tmp_c_mkn = &(tmp_c[addr0]);
				const double *restrict tmp_c_mlkn = tmp_c_mkn + len;
				const double *restrict een_rescaled_n_mnw =
					&(een_rescaled_n[addr1]);
				const double *restrict een_rescaled_n_mlnw =
					een_rescaled_n_mnw + len;
				const double *restrict dtmp_c_mknw = &(dtmp_c[addr0 << 2]);
				const double *restrict dtmp_c_mlknw = dtmp_c_mknw + len4;
				const double *restrict een_rescaled_n_deriv_e_mnw =
					&(een_rescaled_n_deriv_e[addr1 << 2]);
				const double *restrict een_rescaled_n_deriv_e_mlnw =
					een_rescaled_n_deriv_e_mnw + len4;
				for (size_t a = 0; a < (size_t)nucl_num; a++) {
					double cn = c_vector_full[a + n * nucl_num];
					if (cn == 0.0)
						continue;

					const size_t ishift = elec_num * a;
					const size_t ishift4 = ishift << 2;

					const double *restrict tmp_c_amlkn = tmp_c_mlkn + ishift;
					const double *restrict tmp_c_amkn = tmp_c_mkn + ishift;
					const double *restrict een_rescaled_n_amnw =
						een_rescaled_n_mnw + ishift;
					const double *restrict een_rescaled_n_amlnw =
						een_rescaled_n_mlnw + ishift;
					const double *restrict dtmp_c_0amknw =
						dtmp_c_mknw + ishift4;
					const double *restrict dtmp_c_0amlknw =
						dtmp_c_mlknw + ishift4;
					const double *restrict een_rescaled_n_deriv_e_0amnw =
						een_rescaled_n_deriv_e_mnw + ishift4;
					const double *restrict een_rescaled_n_deriv_e_0amlnw =
						een_rescaled_n_deriv_e_mlnw + ishift4;

					const double *restrict dtmp_c_1amknw =
						dtmp_c_0amknw + elec_num;
					const double *restrict dtmp_c_1amlknw =
						dtmp_c_0amlknw + elec_num;
					const double *restrict dtmp_c_2amknw =
						dtmp_c_0amknw + elec_num2;
					const double *restrict dtmp_c_2amlknw =
						dtmp_c_0amlknw + elec_num2;
					const double *restrict dtmp_c_3amknw =
						dtmp_c_0amknw + elec_num3;
					const double *restrict dtmp_c_3amlknw =
						dtmp_c_0amlknw + elec_num3;
					const double *restrict een_rescaled_n_deriv_e_1amnw =
						een_rescaled_n_deriv_e_0amnw + elec_num;
					const double *restrict een_rescaled_n_deriv_e_1amlnw =
						een_rescaled_n_deriv_e_0amlnw + elec_num;
					const double *restrict een_rescaled_n_deriv_e_2amnw =
						een_rescaled_n_deriv_e_0amnw + elec_num2;
					const double *restrict een_rescaled_n_deriv_e_2amlnw =
						een_rescaled_n_deriv_e_0amlnw + elec_num2;
					const double *restrict een_rescaled_n_deriv_e_3amnw =
						een_rescaled_n_deriv_e_0amnw + elec_num3;
					const double *restrict een_rescaled_n_deriv_e_3amlnw =
						een_rescaled_n_deriv_e_0amlnw + elec_num3;
					double *const restrict factor_een_deriv_e_1nw =
						factor_een_deriv_e_0nw + elec_num;
					double *const restrict factor_een_deriv_e_2nw =
						factor_een_deriv_e_0nw + elec_num2;
					double *const restrict factor_een_deriv_e_3nw =
						factor_een_deriv_e_0nw + elec_num3;

					for (size_t j = 0; j < (size_t)elec_num; ++j) {
						factor_een_deriv_e_0nw[j] +=
							cn *
							(tmp_c_amkn[j] * een_rescaled_n_deriv_e_0amlnw[j] +
							 dtmp_c_0amknw[j] * een_rescaled_n_amlnw[j] +
							 dtmp_c_0amlknw[j] * een_rescaled_n_amnw[j] +
							 tmp_c_amlkn[j] * een_rescaled_n_deriv_e_0amnw[j]);
						tmp3[j] =
							dtmp_c_0amknw[j] *
								een_rescaled_n_deriv_e_0amlnw[j] +
							dtmp_c_0amlknw[j] * een_rescaled_n_deriv_e_0amnw[j];
					}

					for (size_t j = 0; j < (size_t)elec_num; ++j) {
						factor_een_deriv_e_1nw[j] +=
							cn *
							(tmp_c_amkn[j] * een_rescaled_n_deriv_e_1amlnw[j] +
							 dtmp_c_1amknw[j] * een_rescaled_n_amlnw[j] +
							 dtmp_c_1amlknw[j] * een_rescaled_n_amnw[j] +
							 tmp_c_amlkn[j] * een_rescaled_n_deriv_e_1amnw[j]);
						tmp3[j] +=
							dtmp_c_1amknw[j] *
								een_rescaled_n_deriv_e_1amlnw[j] +
							dtmp_c_1amlknw[j] * een_rescaled_n_deriv_e_1amnw[j];
					}

					for (size_t j = 0; j < (size_t)elec_num; ++j) {
						factor_een_deriv_e_2nw[j] +=
							cn *
							(tmp_c_amkn[j] * een_rescaled_n_deriv_e_2amlnw[j] +
							 dtmp_c_2amknw[j] * een_rescaled_n_amlnw[j] +
							 dtmp_c_2amlknw[j] * een_rescaled_n_amnw[j] +
							 tmp_c_amlkn[j] * een_rescaled_n_deriv_e_2amnw[j]);
						tmp3[j] +=
							dtmp_c_2amknw[j] *
								een_rescaled_n_deriv_e_2amlnw[j] +
							dtmp_c_2amlknw[j] * een_rescaled_n_deriv_e_2amnw[j];
					}

					for (size_t j = 0; j < (size_t)elec_num; ++j) {
						factor_een_deriv_e_3nw[j] +=
							cn *
							(tmp_c_amkn[j] * een_rescaled_n_deriv_e_3amlnw[j] +
							 dtmp_c_3amknw[j] * een_rescaled_n_amlnw[j] +
							 dtmp_c_3amlknw[j] * een_rescaled_n_amnw[j] +
							 tmp_c_amlkn[j] * een_rescaled_n_deriv_e_3amnw[j] +
							 tmp3[j] * 2.0);
					}
				}
			}
		}
	}
	qmckl_free_device(context, tmp3);
	return info;
}

// Electron/electron/nucleus deriv
qmckl_exit_code_device
qmckl_compute_jastrow_factor_een_rescaled_e_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t cord_num,
	const double rescale_factor_ee, const double *coord_ee,
	const double *ee_distance, const double *een_rescaled_e,
	double *const een_rescaled_e_deriv_e) {

	double x, rij_inv, kappa_l;
	int i, j, k, l, nw, ii;

	double *elec_dist_deriv_e =
		qmckl_malloc_device(context, 4 * elec_num * elec_num * sizeof(double));

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (cord_num < 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

#pragma acc kernels deviceptr(coord_ee, ee_distance, een_rescaled_e,           \
								  een_rescaled_e_deriv_e)
	{
		// Prepare table of exponentiated distances raised to appropriate power
		for (int i = 0; i < elec_num * 4 * elec_num * (cord_num + 1) * walk_num;
			 i++) {
			een_rescaled_e_deriv_e[i] = 0.0;
		}
		for (nw = 0; nw < walk_num; nw++) {
			for (int j = 0; j < elec_num; j++) {
				for (int i = 0; i < elec_num; i++) {
					rij_inv = 1.0 / ee_distance[i + j * elec_num +
												nw * elec_num * elec_num];
					for (int ii = 0; ii < 3; ii++) {
						elec_dist_deriv_e[ii + i * 4 + j * 4 * elec_num] =
							(coord_ee[i + ii * elec_num + nw * elec_num * 3] -
							 coord_ee[j + ii * elec_num + nw * elec_num * 3]) *
							rij_inv;
					}
					elec_dist_deriv_e[3 + i * 4 + j * 4 * elec_num] =
						2.0 * rij_inv;
				}

				for (int ii = 0; ii < 4; ii++) {
					elec_dist_deriv_e[ii + j * 4 + j * 4 * elec_num] = 0.0;
				}
			}

			// prepare the actual een table
			for (l = 0; l < cord_num; l++) {
				kappa_l = -l * rescale_factor_ee;
				for (j = 0; j < elec_num; j++) {
					for (i = 0; i < elec_num; i++) {
						een_rescaled_e_deriv_e[i + 0 * elec_num +
											   j * elec_num * 4 + l * elec_num +
											   nw * elec_num] =
							kappa_l *
							elec_dist_deriv_e[0 + i * 4 + j * 4 * elec_num];
						een_rescaled_e_deriv_e[i + 1 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							kappa_l *
							elec_dist_deriv_e[1 + i * 4 + j * 4 * elec_num];
						een_rescaled_e_deriv_e[i + 2 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							kappa_l *
							elec_dist_deriv_e[2 + i * 4 + j * 4 * elec_num];
						een_rescaled_e_deriv_e[i + 3 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							kappa_l *
							elec_dist_deriv_e[3 + i * 4 + j * 4 * elec_num];

						een_rescaled_e_deriv_e[i + 3 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							een_rescaled_e_deriv_e[i + 3 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] +
							een_rescaled_e_deriv_e[i + 0 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
								een_rescaled_e_deriv_e[i + 0 * elec_num +
													   j * elec_num +
													   l * elec_num +
													   nw * elec_num] +
							een_rescaled_e_deriv_e[i + 1 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
								een_rescaled_e_deriv_e[i + 1 * elec_num +
													   j * elec_num +
													   l * elec_num +
													   nw * elec_num] +
							een_rescaled_e_deriv_e[i + 2 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
								een_rescaled_e_deriv_e[i + 2 * elec_num +
													   j * elec_num +
													   l * elec_num +
													   nw * elec_num];

						een_rescaled_e_deriv_e[i + 0 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							een_rescaled_e_deriv_e[i + 0 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
							een_rescaled_e[i + j + l + nw];
						een_rescaled_e_deriv_e[i + 2 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							een_rescaled_e_deriv_e[i + 1 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
							een_rescaled_e[i + j + l + nw];
						een_rescaled_e_deriv_e[i + 2 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							een_rescaled_e_deriv_e[i + 2 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
							een_rescaled_e[i + j + l + nw];
						een_rescaled_e_deriv_e[i + 3 * elec_num + j * elec_num +
											   l * elec_num + nw * elec_num] =
							een_rescaled_e_deriv_e[i + 3 * elec_num +
												   j * elec_num + l * elec_num +
												   nw * elec_num] *
							een_rescaled_e[i + j + l + nw];
					}
				}
			}
		}
	}
	qmckl_free_device(context, elec_dist_deriv_e);
	return info;
}

qmckl_exit_code_device
qmckl_compute_jastrow_factor_een_rescaled_n_deriv_e_device(
	const qmckl_context_device context, const int64_t walk_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const int64_t cord_num,
	const double *rescale_factor_en, const double *coord_ee,
	const double *coord_en, const double *en_distance,
	const double *een_rescaled_n, double *const een_rescaled_n_deriv_e) {

	double *elnuc_dist_deriv_e =
		qmckl_malloc_device(context, 4 * elec_num * nucl_num * sizeof(double));

	double x, ria_inv, kappa_l;
	int i, a, k, l, nw, ii;

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_4_DEVICE;
		return info;
	}

	if (cord_num < 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

	// Prepare table of exponentiated distances raised to appropriate power
	for (int i = 0; i < elec_num * 4 * nucl_num * (cord_num + 1) * walk_num;
		 i++)
		een_rescaled_n_deriv_e[i] = 0.0;

	for (int nw = 0; nw < walk_num; nw++) {

		// Prepare the actual een table
		for (int a = 0; a < nucl_num; a++) {
			for (int i = 0; i < elec_num; i++) {
				ria_inv =
					1.0 /
					en_distance[a + i * nucl_num + nw * nucl_num * elec_num];
				for (int ii = 0; ii < 3; ii++) {
					elnuc_dist_deriv_e[ii + i + a] =
						(coord_ee[i + ii * elec_num + nw * elec_num * 4] -
						 coord_en[a + ii * nucl_num]) *
						ria_inv;
				}
				elnuc_dist_deriv_e[3 + i + a] = 2.0 * ria_inv;
			}
		}

		for (int l = 0; l < cord_num; l++) {
			for (int a = 0; a < (nucl_num + 1); a++) {
				kappa_l = -((double)l) * rescale_factor_en[type_nucl_vector[a]];
				for (int i = 0; i < elec_num; i++) {
					een_rescaled_n_deriv_e[i + 0 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						kappa_l * elnuc_dist_deriv_e[0 + i + a];
					een_rescaled_n_deriv_e[i + 1 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						kappa_l * elnuc_dist_deriv_e[1 + i + a];
					een_rescaled_n_deriv_e[i + 2 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						kappa_l * elnuc_dist_deriv_e[2 + i + a];
					een_rescaled_n_deriv_e[i + 3 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						kappa_l * elnuc_dist_deriv_e[3 + i + a];

					een_rescaled_n_deriv_e[i + 3 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						een_rescaled_n_deriv_e[i + 3 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] +
						een_rescaled_n_deriv_e[i + 0 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
							een_rescaled_n_deriv_e[i + 1 * elec_num +
												   a * elec_num * 4 +
												   l * elec_num * 4 * nucl_num +
												   nw * elec_num * 4 *
													   nucl_num *
													   (cord_num + 1)] +
						een_rescaled_n_deriv_e[i + 1 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
							een_rescaled_n_deriv_e[i + 2 * elec_num +
												   a * elec_num * 4 +
												   l * elec_num * 4 * nucl_num +
												   nw * elec_num * 4 *
													   nucl_num *
													   (cord_num + 1)] +
						een_rescaled_n_deriv_e[i + 2 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
							een_rescaled_n_deriv_e
								[i + 3 * elec_num + a * elec_num * 4 +
								 l * elec_num * 4 * nucl_num +
								 nw * elec_num * 4 * nucl_num * (cord_num + 1)];

					een_rescaled_n_deriv_e[i + 0 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						een_rescaled_n_deriv_e[i + 0 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
						een_rescaled_n[i + a * elec_num +
									   l * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)];
					een_rescaled_n_deriv_e[i + 1 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						een_rescaled_n_deriv_e[i + 1 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
						een_rescaled_n[i + a * elec_num +
									   l * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)];
					een_rescaled_n_deriv_e[i + 2 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						een_rescaled_n_deriv_e[i + 2 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
						een_rescaled_n[i + a * elec_num +
									   l * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)];
					een_rescaled_n_deriv_e[i + 3 * elec_num + a * elec_num * 4 +
										   l * elec_num * 4 * nucl_num +
										   nw * elec_num * 4 * nucl_num *
											   (cord_num + 1)] =
						een_rescaled_n_deriv_e[i + 3 * elec_num +
											   a * elec_num * 4 +
											   l * elec_num * 4 * nucl_num +
											   nw * elec_num * 4 * nucl_num *
												   (cord_num + 1)] *
						een_rescaled_n[i + a * elec_num +
									   l * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)];
				}
			}
		}
	}

	qmckl_free_device(context, elnuc_dist_deriv_e);
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

#pragma acc kernels deviceptr(coord, ee_distance_rescaled)
	{
		for (int k = 0; k < walk_num; k++) {
			// TODO Fix undefined reference
			/*
			info = qmckl_distance_rescaled_device(
				context, 'T', 'T', elec_num, elec_num, coord + (k * elec_num),
				elec_num * walk_num, coord + (k * elec_num),
				elec_num * walk_num,
				ee_distance_rescaled + (k * elec_num * elec_num), elec_num,
				rescale_factor_ee);
				*/
			if (info != QMCKL_SUCCESS_DEVICE) {
				break;
			}
		}
	}
	return info;
}

qmckl_exit_code_device qmckl_compute_ee_distance_rescaled_deriv_e_device(
	const qmckl_context_device context, const int64_t elec_num,
	const double rescale_factor_ee, const int64_t walk_num, const double *coord,
	double *const ee_distance_rescaled_deriv_e) {

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

	for (int k = 0; k < walk_num; k++) {
		info = qmckl_distance_rescaled_deriv_e_device(
			context, 'T', 'T', elec_num, elec_num, coord + (k * elec_num),
			elec_num * walk_num, coord + (k * elec_num), elec_num * walk_num,
			ee_distance_rescaled_deriv_e + (k * 4 * elec_num * elec_num),
			elec_num, rescale_factor_ee);
		if (info != QMCKL_SUCCESS_DEVICE)
			break;
	}

	return info;
}

qmckl_exit_code_device qmckl_compute_en_distance_rescaled_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t nucl_num, const int64_t type_nucl_num,
	int64_t *const type_nucl_vector, const double *rescale_factor_en,
	const int64_t walk_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance_rescaled) {

	int i, k;
	double *coord = qmckl_malloc_device(context, 3 * sizeof(double));

	qmckl_exit_code_device info = QMCKL_SUCCESS_DEVICE;

	if (context == QMCKL_NULL_CONTEXT_DEVICE) {
		info = QMCKL_INVALID_CONTEXT_DEVICE;
		return info;
	}

	if (elec_num <= 0) {
		info = QMCKL_INVALID_ARG_2_DEVICE;
		return info;
	}

	if (nucl_num <= 0) {
		info = QMCKL_INVALID_ARG_3_DEVICE;
		return info;
	}

	if (walk_num <= 0) {
		info = QMCKL_INVALID_ARG_5_DEVICE;
		return info;
	}

#pragma acc kernels deviceptr(coord, nucl_coord)
	{
		for (i = 0; i < nucl_num; i++) {

			coord[0] = nucl_coord[i + 0 * nucl_num];
			coord[1] = nucl_coord[i + 1 * nucl_num];
			coord[2] = nucl_coord[i + 2 * nucl_num];

			for (k = 0; k < walk_num; k++) {
				// TODO Fix undefined reference
				/*
				info = qmckl_distance_rescaled_device(
					context, 'T', 'T', elec_num, 1,
					elec_coord + k * elec_num + elec_num * walk_num,
					elec_num * walk_num, coord, 1,
					en_distance_rescaled + i * elec_num +
						k * elec_num * nucl_num,
					elec_num, rescale_factor_en[type_nucl_vector[i]]);
				*/
				if (info != QMCKL_SUCCESS_DEVICE) {
					break;
				}
			}
		}
	}

	qmckl_free_device(context, coord);
	return info;
}

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

	const size_t elec_pairs = (size_t)(elec_num * (elec_num - 1)) / 2;
	const size_t len_een_ij = (size_t)elec_pairs * (cord_num + 1);
	double *een_rescaled_e_ij =
		qmckl_malloc_device(context, len_een_ij * sizeof(double));
#pragma acc kernels deviceptr(ee_distance, een_rescaled_e, een_rescaled_e_ij)
	{
		// Prepare table of exponentiated distances raised to appropriate power
		// init
		for (int i = 0; i < walk_num * (cord_num + 1) * elec_num * elec_num;
			 i++)
			een_rescaled_e[i] = 0;

		// number of elements for the een_rescaled_e_ij[N_e*(N_e-1)/2][cord+1]
		// probably in C is better [cord+1, Ne*(Ne-1)/2]
		// elec_pairs = (elec_num * (elec_num - 1)) / 2;
		// len_een_ij = elec_pairs * (cord_num + 1);
		const size_t e2 = elec_num * elec_num;

		for (size_t nw = 0; nw < (size_t)walk_num; ++nw) {

			for (size_t kk = 0; kk < len_een_ij; ++kk) {
				een_rescaled_e_ij[kk] = 0.0;
			}
			for (size_t kk = 0; kk < elec_pairs; ++kk) {
				een_rescaled_e_ij[kk] = 1.0;
			}

			size_t kk = 0;
			for (size_t i = 0; i < (size_t)elec_num; ++i) {
				for (size_t j = 0; j < i; ++j) {
					een_rescaled_e_ij[j + kk + elec_pairs] =
						-rescale_factor_ee *
						ee_distance[j + i * elec_num + nw * e2];
				}
				kk += i;
			}

			for (size_t k = elec_pairs; k < 2 * elec_pairs; ++k) {
				een_rescaled_e_ij[k] = exp(een_rescaled_e_ij[k]);
			}

			for (size_t l = 2; l < (size_t)(cord_num + 1); ++l) {
				for (size_t k = 0; k < elec_pairs; ++k) {
					// een_rescaled_e_ij(k, l + 1) = een_rescaled_e_ij(k, l + 1
					// - 1)
					// * een_rescaled_e_ij(k, 2)
					een_rescaled_e_ij[k + l * elec_pairs] =
						een_rescaled_e_ij[k + (l - 1) * elec_pairs] *
						een_rescaled_e_ij[k + elec_pairs];
				}
			}

			double *const een_rescaled_e_ =
				&(een_rescaled_e[nw * (cord_num + 1) * e2]);
			// prepare the actual een table
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
	}
	qmckl_free_device(context, een_rescaled_e_ij);
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

#pragma acc kernels deviceptr(type_nucl_vector, rescale_factor_en, en_distance)
	{
		// Prepare table of exponentiated distances raised to appropriate power
		for (int i = 0; i < (walk_num * (cord_num + 1) * nucl_num * elec_num);
			 ++i) {
			een_rescaled_n[i] = 1.0;
		}

		for (int nw = 0; nw < walk_num; ++nw) {
			for (int a = 0; a < nucl_num; ++a) {
				for (int i = 0; i < elec_num; ++i) {
					een_rescaled_n[i + a * elec_num +
								   nw * elec_num * nucl_num * (cord_num + 1)] =
						1.0;
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
						een_rescaled_n[i + a * elec_num +
									   l * elec_num * nucl_num +
									   nw * elec_num * nucl_num *
										   (cord_num + 1)] =
							een_rescaled_n[i + a * elec_num +
										   (l - 1) * elec_num * nucl_num +
										   nw * elec_num * nucl_num *
											   (cord_num + 1)] *
							een_rescaled_n[i + a * elec_num +
										   elec_num * nucl_num +
										   nw * elec_num * nucl_num *
											   (cord_num + 1)];
					}
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

#pragma acc kernels deviceptr(type_nucl_vector, c_vector, c_vector_full)
	{
		for (int i = 0; i < dim_c_vector; ++i) {
			for (int a = 0; a < nucl_num; ++a) {
				c_vector_full[a + i * nucl_num] =
					c_vector[(type_nucl_vector[a] - 1) + i * type_nucl_num];
			}
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

#pragma acc kernels deviceptr(lkpm_combined_index)
	{
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

#pragma acc kernels deviceptr(een_rescaled_e, een_rescaled_n)
	{
		for (int64_t nw = 0; nw < walk_num; ++nw) {
			for (int64_t i = 0; i < cord_num; ++i) {
				/* TODO Replace by a device BLAS call or write a temporary
				manual DGEMM info = qmckl_dgemm(context, TransA, TransB, M, N,
				K, alpha,
								&(een_rescaled_e[af * (i + nw * (cord_num +
				1))]), LDA, &(een_rescaled_n[bf * nw]), LDB, beta,
								&(tmp_c[cf * (i + nw * cord_num)]), LDC);
				 */
			}
		}
	}
	return info;
}

qmckl_exit_code_device qmckl_compute_dtmp_c_device(
	const qmckl_context_device context, const int64_t cord_num,
	const int64_t elec_num, const int64_t nucl_num, const int64_t walk_num,
	const double *een_rescaled_e_deriv_e, const double *een_rescaled_n,
	double *const dtmp_c) {

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

	const int64_t M = 4 * elec_num;
	const int64_t N = nucl_num * (cord_num + 1);
	const int64_t K = elec_num;

	const int64_t LDA = 4 * elec_num;
	const int64_t LDB = elec_num;
	const int64_t LDC = 4 * elec_num;

	const int64_t af = elec_num * elec_num * 4;
	const int64_t bf = elec_num * nucl_num * (cord_num + 1);
	const int64_t cf = elec_num * 4 * nucl_num * (cord_num + 1);

#pragma acc kernels deviceptr(een_rescaled_e_deriv_e, een_rescaled_n, dtmp_c)
	{
		for (int64_t nw = 0; nw < walk_num; ++nw) {
			for (int64_t i = 0; i < cord_num; ++i) {
				/* TODO Replace by a device BLAS call or write a temporary
				  manual DGEMM info = qmckl_dgemm(context, TransA, TransB, M, N,
				  K, alpha,
									 &(een_rescaled_e_deriv_e[af*(i+nw*(cord_num+1))]),
				  LDA,
									 &(een_rescaled_n[bf*nw]), LDB, beta,
									 &(dtmp_c[cf*(i+nw*cord_num)]), LDC);
							   */
			}
		}
	}

	return info;
}
