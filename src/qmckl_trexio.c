#include "../include/qmckl_trexio.h"

//**********
// TREXIO CONTEXT FILL
//**********

qmckl_exit_code_device
qmckl_trexio_read_electron_X_device(qmckl_context_device context,
									trexio_t *file) {

	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	int rcio = 0;

	int64_t up_num = 0L;
	int64_t dn_num = 0L;

	rcio = trexio_read_electron_up_num_64(file, &up_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_electron_up_num",
									 trexio_string_of_error(rcio));
	}

	assert(up_num >= 0L);

	rcio = trexio_read_electron_dn_num_64(file, &dn_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_electron_dn_num",
									 trexio_string_of_error(rcio));
	}

	assert(dn_num >= 0L);

	qmckl_exit_code_device rc;
	rc = qmckl_set_electron_num_device(context, up_num, dn_num);
	return rc;
}

qmckl_exit_code_device
qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
								   trexio_t *file) {
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	qmckl_exit_code_device rc;
	int rcio = 0;

	int64_t nucleus_num = 0L;

	rcio = trexio_read_nucleus_num_64(file, &nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_nucleus_num",
									 trexio_string_of_error(rcio));
	}

	assert(nucleus_num > 0);
	rc = qmckl_set_nucleus_num_device(context, nucleus_num);

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;
		mem_info.size = nucleus_num * sizeof(double);

		double *nucl_charge_h = (double *)qmckl_malloc_host(context, mem_info);
		double *nucl_charge_d = (double *)qmckl_malloc_device(
			context, nucleus_num * sizeof(double));

		if (nucl_charge_h == NULL || nucl_charge_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
				"qmckl_trexio_read_nucleus_X_device", NULL);
		}

		assert(nucl_charge_h != NULL && nucl_charge_d != NULL);

		rcio = trexio_read_safe_nucleus_charge_64(file, nucl_charge_h,
												  nucleus_num);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	qmckl_memory_info_struct_device mem_info =
		qmckl_memory_info_struct_zero_device;
	mem_info.size = nucleus_num * 3 * sizeof(double);

	double *nucl_coord_h = (double *)qmckl_malloc_host(context, mem_info);
	double *nucl_coord_d = (double *)qmckl_malloc_device(
		context, nucleus_num * 3 * sizeof(double));

	if (nucl_coord_h == NULL || nucl_coord_d == NULL) {
		return qmckl_failwith_device(context, QMCKL_ALLOCATION_FAILED_DEVICE,
									 "qmckl_trexio_read_nucleus_X_device",
									 NULL);
	}

	assert(nucl_coord_h != NULL && nucl_coord_d != NULL);

	rcio =
		trexio_read_safe_nucleus_coord_64(file, nucl_coord_h, 3 * nucleus_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_nucleus_charge",
									 trexio_string_of_error(rcio));
	}

	qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord_h, mem_info.size);
	rc = qmckl_set_nucleus_coord_device(context, 'N', nucl_coord_d,
										3 * nucleus_num);

	qmckl_free_host(context, nucl_coord_h);
	qmckl_free_device(context, nucl_coord_d);
	nucl_coord_h = NULL;
	nucl_coord_d = NULL;

	if (rc != QMCKL_SUCCESS_DEVICE) {
		return rc;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_trexio_read_ao_X_device(qmckl_context_device context, trexio_t *file) {
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	qmckl_exit_code_device rc;
	int rcio = 0;
	int64_t nucleus_num = 0L;

	rc = qmckl_get_nucleus_num_device(context, &nucleus_num);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

#define MAX_STR_LEN 1024
	char basis_type[MAX_STR_LEN];

	rcio = trexio_read_basis_type(file, basis_type, MAX_STR_LEN);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_basis_type",
									 trexio_string_of_error(rcio));
	}

	if (basis_type[0] == 'G') {
		rc = qmckl_set_ao_basis_type_device(context, basis_type[0]);
	} else {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_basis_type",
									 "Invalid basis type");
	}

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t shell_num = 0L;

	rcio = trexio_read_basis_shell_num_64(file, &shell_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_basis_shell_num",
									 trexio_string_of_error(rcio));
	}

	assert(shell_num > 0);
	rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t prim_num = 0L;

	rcio = trexio_read_basis_prim_num_64(file, &prim_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_basis_prim_num",
									 trexio_string_of_error(rcio));
	}

	assert(prim_num > 0);
	rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t ao_num = 0LL;

	rcio = trexio_read_ao_num_64(file, &ao_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_ao_num",
									 trexio_string_of_error(rcio));
	}

	assert(ao_num > 0);
	rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;
		int64_t *nucleus_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_index_d = (int64_t *)qmckl_malloc_device(
			context, nucleus_num * sizeof(int64_t));

		if (nucleus_index_h == NULL || nucleus_index_d == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
										 "trexio_read_basis_nucleus_index",
										 trexio_string_of_error(rcio));
		}

		/* Reformat data */
		rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d,
													 nucleus_num);
		if (rc != QMCKL_SUCCESS_DEVICE) {
			qmckl_free_host(context, nucleus_index_h);
			qmckl_free_device(context, nucleus_index_d);
			nucleus_index_h = NULL;
			nucleus_index_d = NULL;
			return rc;
		}

		for (int i = shell_num - 1; i >= 0; --i) {
			int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_index_h);
				qmckl_free_device(context, nucleus_index_d);
				nucleus_index_h = NULL;
				nucleus_index_d = NULL;
				return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = nucleus_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *nucleus_shell_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *nucleus_shell_num_d = (int64_t *)qmckl_malloc_device(
			context, nucleus_num * sizeof(int64_t));

		if (nucleus_shell_num_h == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
										 "trexio_read_basis_nucleus_shell_num",
										 trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < nucleus_num; ++i) {
			nucleus_shell_num_h[i] = 0;
		}

		for (int i = 0; i < shell_num; ++i) {
			int k = tmp_array[i];
			if (k < 0 || k >= nucleus_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, nucleus_shell_num_h);
				qmckl_free_device(context, nucleus_shell_num_d);
				nucleus_shell_num_h = NULL;
				return qmckl_failwith_device(
					context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int32_t);

		int32_t *shell_ang_mom_h =
			(int32_t *)qmckl_malloc_host(context, mem_info);
		int32_t *shell_ang_mom_d = (int32_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int32_t));

		if (shell_ang_mom_h == NULL || shell_ang_mom_d == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_num_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_num_d = (int64_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int64_t));

		if (shell_prim_num_h == NULL || shell_prim_num_d == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
										 "trexio_read_basis_shell_prim_num",
										 trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = 0; i < shell_num; ++i) {
			shell_prim_num_h[i] = 0;
		}

		for (int i = 0; i < prim_num; ++i) {
			int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				qmckl_free_host(context, shell_prim_num_h);
				qmckl_free_device(context, shell_prim_num_d);
				return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(int64_t);
		size_t size_backup = mem_info.size;

		int64_t *shell_prim_index_h =
			(int64_t *)qmckl_malloc_host(context, mem_info);
		int64_t *shell_prim_index_d = (int64_t *)qmckl_malloc_device(
			context, shell_num * sizeof(int64_t));

		if (shell_prim_index_h == NULL || shell_prim_index_d == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
										 "qmckl_trexio_read_basis_shell_"
										 "prim_index_X_device",
										 NULL);
		}

		assert(shell_prim_index_h != NULL && shell_prim_index_d != NULL);

		/* Allocate temporary array */
		mem_info.size = prim_num * sizeof(int64_t);

		int64_t *tmp_array = (int64_t *)qmckl_malloc_host(context, mem_info);

		if (tmp_array == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
										 "trexio_read_basis_shell_prim_index",
										 trexio_string_of_error(rcio));
		}

		/* Reformat data */
		for (int i = prim_num - 1; i >= 0; --i) {
			int k = tmp_array[i];
			if (k < 0 || k >= shell_num) {
				qmckl_free_host(context, tmp_array);
				tmp_array = NULL;
				qmckl_free_host(context, shell_prim_index_h);
				qmckl_free_device(context, shell_prim_index_d);
				shell_prim_index_h = NULL;
				shell_prim_index_d = NULL;
				return qmckl_failwith_device(
					context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = shell_num * sizeof(double);

		double *shell_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *shell_factor_d =
			(double *)qmckl_malloc_device(context, shell_num * sizeof(double));

		if (shell_factor_h == NULL || shell_factor_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *exponent_h = (double *)qmckl_malloc_host(context, mem_info);
		double *exponent_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (exponent_h == NULL || exponent_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *coefficient_h = (double *)qmckl_malloc_host(context, mem_info);
		double *coefficient_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (coefficient_h == NULL || coefficient_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = prim_num * sizeof(double);

		double *prim_factor_h = (double *)qmckl_malloc_host(context, mem_info);
		double *prim_factor_d =
			(double *)qmckl_malloc_device(context, prim_num * sizeof(double));

		if (prim_factor_h == NULL || prim_factor_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;

		/* Allocate array for data */
		mem_info.size = ao_num * sizeof(double);

		double *ao_normalization_h =
			(double *)qmckl_malloc_host(context, mem_info);
		double *ao_normalization_d =
			(double *)qmckl_malloc_device(context, ao_num * sizeof(double));

		if (ao_normalization_h == NULL || ao_normalization_d == NULL) {
			return qmckl_failwith_device(
				context, QMCKL_ALLOCATION_FAILED_DEVICE,
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
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
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

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	return QMCKL_SUCCESS_DEVICE;
}

qmckl_exit_code_device
qmckl_trexio_read_mo_X_device(qmckl_context_device context, trexio_t *file) {
	assert(context != (qmckl_context_device)0);
	assert(file != NULL);

	qmckl_exit_code_device rc;
	int rcio = 0;
	int64_t ao_num = 0L;

	rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num);
	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	int64_t mo_num = 0L;

	rcio = trexio_read_mo_num_64(file, &mo_num);
	if (rcio != TREXIO_SUCCESS) {
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "trexio_read_mo_num",
									 trexio_string_of_error(rcio));
	}

	assert(mo_num > 0);
	rc = qmckl_set_mo_basis_mo_num_device(context, mo_num);

	if (rc != QMCKL_SUCCESS_DEVICE)
		return rc;

	{
		qmckl_memory_info_struct_device mem_info =
			qmckl_memory_info_struct_zero_device;
		mem_info.size = ao_num * mo_num * sizeof(double);

		double *mo_coef_h = (double *)qmckl_malloc_host(context, mem_info);
		double *mo_coef_d = (double *)qmckl_malloc_device(
			context, ao_num * mo_num * sizeof(double));

		if (mo_coef_h == NULL || mo_coef_d == NULL) {
			return qmckl_failwith_device(context,
										 QMCKL_ALLOCATION_FAILED_DEVICE,
										 "qmckl_trexio_read_mo_X_device", NULL);
		}

		assert(mo_coef_h != NULL && mo_coef_d != NULL);

		rcio = trexio_read_mo_coefficient_64(file, mo_coef_h);
		if (rcio != TREXIO_SUCCESS) {
			return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
										 "trexio_read_mo_coefficient",
										 trexio_string_of_error(rcio));
		}

		qmckl_memcpy_H2D(context, mo_coef_d, mo_coef_h, mem_info.size);
		rc = qmckl_set_mo_basis_coefficient_device(context, mo_coef_d);

		qmckl_free_host(context, mo_coef_h);
		qmckl_free_device(context, mo_coef_d);
		mo_coef_h = NULL;
		mo_coef_d = NULL;

		if (rc != QMCKL_SUCCESS_DEVICE)
			return rc;
	}

	return QMCKL_SUCCESS_DEVICE;
}

trexio_t *qmckl_trexio_open_X_device(char *file_name,
									 qmckl_exit_code_device *rc) {
	*rc = QMCKL_SUCCESS_DEVICE;
	trexio_t *file = NULL;

	file = trexio_open(file_name, 'r', TREXIO_TEXT, rc);
	if (file != NULL)
		return file;

	file = trexio_open(file_name, 'r', TREXIO_HDF5, rc);
	if (file != NULL)
		return file;

	*rc = QMCKL_FAILURE_DEVICE;
	/* TODO
	  return qmckl_failwith( context,
							 QMCKL_FAILURE,
							 "trexio_read_electron_up_num",
							 trexio_string_of_error(rcio));
							 */
	return NULL;
}


qmckl_exit_code_device qmckl_trexio_read_device(qmckl_context_device context,
												char *file_name,
												int64_t size_max) {
	if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_exit_code_device rc;
	char file_name_new[size_max + 1];
	strncpy(file_name_new, file_name, size_max + 1);
	file_name_new[size_max] = '\0';

	trexio_t *file = qmckl_trexio_open_X_device(file_name_new, &rc);
	if (file == NULL) {
		trexio_close(file);
		return qmckl_failwith_device(context, QMCKL_INVALID_ARG_2_DEVICE,
									 "qmckl_trexio_read_device",
									 trexio_string_of_error(rc));
	}

	assert(file != NULL);

	rc = qmckl_trexio_read_electron_X_device(context, file);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		trexio_close(file);
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_trexio_read_device",
									 "Error reading electron");
	}

	rc = qmckl_trexio_read_nucleus_X_device(context, file);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		trexio_close(file);
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_trexio_read_device",
									 "Error reading nucleus");
	}

	rc = qmckl_trexio_read_ao_X_device(context, file);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		trexio_close(file);
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_trexio_read_device",
									 "Error reading AOs");
	}

	rc = qmckl_trexio_read_mo_X_device(context, file);
	if (rc != QMCKL_SUCCESS_DEVICE) {
		trexio_close(file);
		return qmckl_failwith_device(context, QMCKL_FAILURE_DEVICE,
									 "qmckl_trexio_omp_read",
									 "Error reading MOs");
	}

	trexio_close(file);
	file = NULL;
	return rc;
}
