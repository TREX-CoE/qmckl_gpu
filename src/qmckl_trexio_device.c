#include "../include/qmckl_trexio_device.h"

//**********
// ADDITIONAL GETTERS/SETTERS
//**********


qmckl_exit_code
qmckl_get_nucleus_num_device (const qmckl_context_device context, int64_t* const num) {

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (num == NULL) {
    return qmckl_failwith_device( context,
                                   QMCKL_INVALID_ARG_2,
                                   "qmckl_get_nucleus_num_device",
                                   "num is a null pointer");
  }

  qmckl_context_struct_device* const ctx = (qmckl_context_struct_device*) context;

  if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_nucleus_num_device",
                             NULL);
  }

  assert (ctx != NULL);

  assert (ctx->nucleus.num >= (int64_t) 0);
  *num = ctx->nucleus.num;

  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_nucleus_num_device(qmckl_context_device context,
                             const int64_t num)
{

  int32_t mask = 1 << 0;

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith_device( context,
                                    QMCKL_NULL_CONTEXT,
                                    "qmckl_set_nucleus_num_device",
                                    NULL);
  }

  qmckl_context_struct_device* const ctx = (qmckl_context_struct_device*) context;


  if (num <= 0) {
    return qmckl_failwith_device( context,
                                  QMCKL_INVALID_ARG_2,
                                  "qmckl_set_nucleus_num",
                                  "num <= 0");
  }

  ctx->nucleus.num = num;

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

  return QMCKL_SUCCESS;

}


qmckl_exit_code
qmckl_set_nucleus_charge_device(qmckl_context_device context,
                                const double* charge,
                                const int64_t size_max,
                                int device_id)
{

  int32_t mask = 1 << 1;

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith_device ( context,
                                     QMCKL_NULL_CONTEXT,
                                     "qmckl_set_nucleus_charge_device",
                                     NULL);
  }


  // This accepts a host pointer and copies it in context as a device pointer

  if (charge == NULL) {
    return qmckl_failwith_device( context,
                                  QMCKL_INVALID_ARG_2,
                                  "qmckl_set_nucleus_charge_device",
                                  "charge is a null pointer");
  }

  int64_t num;
  qmckl_exit_code rc;

  rc = qmckl_get_nucleus_num_device(context, &num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (num > size_max) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_nucleus_charge_device",
                           "Array too small");
  }

  ctx->nucleus.charge = qmckl_vector_alloc_device(context, num, device_id);
  rc = qmckl_vector_of_double_device(context, charge, num, &(ctx->nucleus.charge), device_id);

  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_nucleus_charge_device",
                           "Error in vector->double* conversion");
  }

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_nucleus_coord_device (qmckl_context_device context,
                                const char transp,
                                const double* coord,
                                const int64_t size_max,
                                int32_t device_id)
{
  int32_t mask = 1 << 2;

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith_device ( context,
                                     QMCKL_NULL_CONTEXT,
                                     "qmckl_set_nucleus_coord_device",
                                     NULL);
  }

  qmckl_exit_code rc;

  const int64_t nucl_num = (int64_t) ctx->nucleus.num;

  if (ctx->nucleus.coord.data_device != NULL) {
    rc = qmckl_matrix_free_device(context, &(ctx->nucleus.coord), device_id);
    if (rc != QMCKL_SUCCESS) return rc;
  }

  ctx->nucleus.coord = qmckl_matrix_alloc_device(context, nucl_num, 3, device_id);

  #pragma use_device_ptr(ctx->nucleus.coord.data)
  {
  if (ctx->nucleus.coord.data_device == NULL) {
    return qmckl_failwith_device( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_nucleus_coord_device",
                           NULL);
  }
  }

  if (size_max < 3*nucl_num) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_set_nucleus_coord_device",
                           "Array too small");
  }

  if (transp == 'N') {
    qmckl_matrix At;
    At = qmckl_matrix_alloc_device(context, 3, nucl_num, device_id);
    rc = qmckl_matrix_of_double_device(context, coord, 3*nucl_num, &At, device_id);
    if (rc != QMCKL_SUCCESS) return rc;
    rc = qmckl_transpose_device(context, At, ctx->nucleus.coord);
  } else {
    rc = qmckl_matrix_of_double_device(context, coord, nucl_num*3, &(ctx->nucleus.coord), device_id);
  }
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);

  return QMCKL_SUCCESS;
}



qmckl_exit_code
qmckl_set_ao_basis_type_device(qmckl_context_device context,
                        const char basis_type)
{
  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_set_ao_basis_type_device",
                           NULL);
   }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (basis_type != 'G' && basis_type != 'S') {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_type_device",
                           NULL);
  }

  int32_t mask = 1;

  ctx->ao_basis.type = basis_type;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis_device(context, device_id);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }

  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_ao_basis_shell_num_device (qmckl_context_device context,
                              const int64_t shell_num)
{
  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_set_ao_basis_shell_num_device",
                           NULL);
   }

    if (shell_num <= 0) {
      return qmckl_failwith_device( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_set_ao_basis_shell_num_device",
                             "shell_num <= 0");
    }

  const int64_t prim_num = ctx->ao_basis.prim_num;

  if (0L < prim_num && prim_num < shell_num) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num_device",
                           "shell_num > prim_num");
  }

  int32_t mask = 1 << 1;

  ctx->ao_basis.shell_num = shell_num;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis_device(context, device_id);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }

  return QMCKL_SUCCESS;

}


// NOTE HERE
qmckl_exit_code
qmckl_set_ao_basis_prim_num (qmckl_context_device context,
                             const int64_t prim_num)
{

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_set_ao_basis_prim_num_device",
                           NULL);
   }


    if (prim_num <= 0) {
      return qmckl_failwith_device( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_set_ao_basis_shell_num_device",
                             "prim_num must be positive");
    }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num <= 0L) {
    return qmckl_failwith_device( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_num_device",
                           "shell_num is not set");
  }

  if (prim_num < shell_num) {
    return qmckl_failwith_device( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num_device",
                           "prim_num < shell_num");
  }

  int32_t mask = 1 << 2;

  ctx->ao_basis.prim_num = prim_num;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis_device(context, device_id);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }

  return QMCKL_SUCCESS;

}



//**********
// CONTEXT FILL
//**********


qmckl_exit_code
qmckl_trexio_read_nucleus_X_device(qmckl_context_device context, trexio_t* const file, int device_id)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;

int64_t nucleus_num = 0L;

rcio = trexio_read_nucleus_num_64(file, &nucleus_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith_device( context,
                         QMCKL_FAILURE,
                         "trexio_read_nucleus_num",
                         trexio_string_of_error(rcio));
}

assert (nucleus_num > 0);
rc = qmckl_set_nucleus_num_device(context, nucleus_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucleus_num * sizeof(double);

  double* nucl_charge = (double*) qmckl_malloc_host(context, mem_info);

  if (nucl_charge == NULL) {
    return qmckl_failwith_device( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_nucleus_X_device",
                           NULL);
  }

  assert (nucl_charge != NULL);

  rcio = trexio_read_safe_nucleus_charge_64(file, nucl_charge, nucleus_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_nucleus_charge",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_nucleus_charge_device(context, nucl_charge, nucleus_num, device_id);

  qmckl_free_host(context, nucl_charge);
  nucl_charge = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;

}

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucleus_num * 3 * sizeof(double);

  double* nucl_coord = (double*) qmckl_malloc_host(context, mem_info);

  if (nucl_coord == NULL) {
    return qmckl_failwith_device( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_nucleus_X_device",
                           NULL);
  }

  assert (nucl_coord != NULL);

  rcio = trexio_read_safe_nucleus_coord_64(file, nucl_coord, 3*nucleus_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith_device( context,
                           QMCKL_FAILURE,
                           "trexio_read_nucleus_charge",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_nucleus_coord_device(context, 'N', nucl_coord, 3*nucleus_num, device_id);

  qmckl_free_host(context, nucl_coord);
  nucl_coord = NULL;

  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  return QMCKL_SUCCESS;
}




qmckl_exit_code
qmckl_trexio_read_ao_X_device(qmckl_context context, trexio_t* const file, int device_id)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;
  int64_t nucleus_num = 0L;

  rc = qmckl_get_nucleus_num_device(context, &nucleus_num);
  if (rc != QMCKL_SUCCESS)
    return rc;

#define MAX_STR_LEN 1024
  char basis_type[MAX_STR_LEN];

  rcio = trexio_read_basis_type(file, basis_type, MAX_STR_LEN);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_type",
                           trexio_string_of_error(rcio));
  }

  if (basis_type[0] == 'G') {
    rc = qmckl_set_ao_basis_type_device(context, basis_type[0]);
  } else {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_type",
                           "Invalid basis type");
  }

  if (rc != QMCKL_SUCCESS)
    return rc;

int64_t shell_num = 0L;

rcio = trexio_read_basis_shell_num_64(file, &shell_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_basis_shell_num",
                         trexio_string_of_error(rcio));
}

assert (shell_num > 0);
rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);

if (rc != QMCKL_SUCCESS)
  return rc;

int64_t prim_num = 0L;

rcio = trexio_read_basis_prim_num_64(file, &prim_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_basis_prim_num",
                         trexio_string_of_error(rcio));
}

assert (prim_num > 0);
// NOTE HERE
rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);

if (rc != QMCKL_SUCCESS)
  return rc;

int64_t ao_num = 0LL;

rcio = trexio_read_ao_num_64(file, &ao_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_ao_num",
                         trexio_string_of_error(rcio));
}

assert (ao_num > 0);
// TODO Device
rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = nucleus_num * sizeof(int64_t);
  int64_t* nucleus_index = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (nucleus_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_index_X_device",
                           NULL);
  }

  assert (nucleus_index != NULL);

  /* Allocate temporary array */
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* tmp_array = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free_host(context, nucleus_index);
    nucleus_index = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_index_X_device",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read in the temporary array */
  rcio = trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, tmp_array);
    tmp_array = NULL;
    qmckl_free_host(context, nucleus_index);
    nucleus_index = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_nucleus_index",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index, nucleus_num, device_id);
  if (rc != QMCKL_SUCCESS) {
    qmckl_free_host(context, nucleus_index);
    nucleus_index = NULL;
    return rc;
  }

  for (int i=shell_num-1 ; i>=0 ; --i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= nucleus_num) {
      qmckl_free_host(context, tmp_array);
      tmp_array = NULL;
      qmckl_free_host(context, nucleus_index);
      nucleus_index = NULL;
      return qmckl_failwith( context,
                              QMCKL_FAILURE,
                              "trexio_read_basis_nucleus_index",
                              "Irrelevant data in TREXIO file");
    }
    nucleus_index[k] = i;
  }

  qmckl_free_host(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index, shell_num, device_id);

  qmckl_free_host(context, nucleus_index);
  nucleus_index = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = nucleus_num * sizeof(int64_t);
  int64_t* nucleus_shell_num = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (nucleus_shell_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_shell_num_X_device",
                           NULL);
  }

  assert (nucleus_shell_num != NULL);

  /* Allocate temporary array */
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* tmp_array = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free_host(context, nucleus_shell_num);
    nucleus_shell_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_shell_num_X",
                           NULL);
  }

  assert (tmp_array != NULL);


  /* Read in the temporary array */
  rcio = trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, tmp_array);
    tmp_array = NULL;
    qmckl_free_host(context, nucleus_shell_num);
    nucleus_shell_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_nucleus_shell_num",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=0 ; i<nucleus_num ; ++i) {
    nucleus_shell_num[i] = 0;
  }

  for (int i=0 ; i<shell_num ; ++i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= nucleus_num) {
      qmckl_free_host(context, tmp_array);
      tmp_array = NULL;
      qmckl_free_host(context, nucleus_shell_num);
      nucleus_shell_num = NULL;
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "trexio_read_basis_nucleus_shell_num",
                             "Irrelevant data in TREXIO file");
    }
    nucleus_shell_num[k] += 1;
  }

  qmckl_free_host(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_nucleus_shell_num_device(context, nucleus_shell_num, shell_num, device_id);

  qmckl_free_host(context, nucleus_shell_num);
  nucleus_shell_num = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int32_t);

  int32_t* shell_ang_mom = (int32_t*) qmckl_malloc_host(context, mem_info);

  if (shell_ang_mom == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_ang_mom_X",
                           NULL);
  }

  assert (shell_ang_mom != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_ang_mom_32(file, shell_ang_mom, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, shell_ang_mom);
    shell_ang_mom = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_ang_mom",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom, shell_num, device_id);

  qmckl_free_host(context, shell_ang_mom);
  shell_ang_mom = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int64_t);

  int64_t* shell_prim_num = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (shell_prim_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_num_X_device",
                           NULL);
  }

  assert (shell_prim_num != NULL);

  /* Allocate temporary array */
  mem_info.size = prim_num * sizeof(int64_t);

  int64_t* tmp_array = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free_host(context, shell_prim_num);
    shell_prim_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_num_X",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_index_64 (file, tmp_array, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, shell_prim_num);
    shell_prim_num = NULL;
    qmckl_free_host(context, tmp_array);
    tmp_array = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_prim_num",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=0 ; i<shell_num ; ++i) {
    shell_prim_num[i] = 0;
  }

  for (int i=0 ; i<prim_num ; ++i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= shell_num) {
      qmckl_free_host(context, tmp_array);
      qmckl_free_host(context, shell_prim_num);
      return qmckl_failwith( context,
                              QMCKL_FAILURE,
                              "trexio_read_basis_shell_prim_num",
                              "Irrelevant data in TREXIO file");
    }
    shell_prim_num[k] += 1;
  }

  qmckl_free_host(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num, shell_num, device_id);

  qmckl_free_host(context, shell_prim_num);
  shell_prim_num = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int64_t);

  int64_t* shell_prim_index = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (shell_prim_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_index_X",
                           NULL);
  }

  assert (shell_prim_index != NULL);

  /* Allocate temporary array */
  mem_info.size = prim_num * sizeof(int64_t);

  int64_t* tmp_array = (int64_t*) qmckl_malloc_host(context, mem_info);

  if (tmp_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_index_X",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, shell_prim_index);
    shell_prim_index = NULL;
    qmckl_free_host(context, tmp_array);
    tmp_array = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_prim_index",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=prim_num-1 ; i>=0 ; --i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= shell_num) {
      qmckl_free_host(context, tmp_array);
      tmp_array = NULL;
      qmckl_free_host(context, shell_prim_index);
      shell_prim_index = NULL;
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "trexio_read_basis_shell_prim_index",
                             "Irrelevant data in TREXIO file");
    }
    shell_prim_index[k] = i;
  }

  qmckl_free_host(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index, shell_num, device_id);

  qmckl_free_host(context, shell_prim_index);
  shell_prim_index = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(double);

  double* shell_factor = (double*) qmckl_malloc_host(context, mem_info);

  if (shell_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_factor_X_device",
                           NULL);
  }

  assert (shell_factor != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_factor_64(file, shell_factor, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, shell_factor);
    shell_factor = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_factor",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor, shell_num, device_id);

  qmckl_free_host(context, shell_factor);
  shell_factor = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* exponent = (double*) qmckl_malloc_host(context, mem_info);

  if (exponent == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_exponent_X",
                           NULL);
  }

  assert (exponent != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_exponent_64(file, exponent, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, exponent);
    exponent = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_exponent",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_exponent_device(context, exponent, prim_num, device_id);

  qmckl_free_host(context, exponent);
  exponent = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* coefficient = (double*) qmckl_malloc_host(context, mem_info);

  if (coefficient == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_coefficient_X_device",
                           NULL);
  }

  assert (coefficient != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_coefficient_64(file, coefficient, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, coefficient);
    coefficient = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_coefficient",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_coefficient_device(context, coefficient, prim_num, device_id);

  qmckl_free_host(context, coefficient);
  coefficient = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* prim_factor = (double*) qmckl_malloc_host(context, mem_info);

  if (prim_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_prim_factor_X_device",
                           NULL);
  }

  assert (prim_factor != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_prim_factor_64(file, prim_factor, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, prim_factor);
    prim_factor = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_prim_factor",
                           trexio_string_of_error(rcio));
  }

  /* Read data */
  rc = qmckl_set_ao_basis_prim_factor_device(context, prim_factor, prim_num, device_id);

  qmckl_free_host(context, prim_factor);
  prim_factor = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = ao_num * sizeof(double);

  double* ao_normalization = (double*) qmckl_malloc_host(context, mem_info);

  if (ao_normalization == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_ao_normalization_X_device",
                           NULL);
  }

  assert (ao_normalization != NULL);

  /* Read data */
  rcio = trexio_read_safe_ao_normalization_64(file, ao_normalization, ao_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free_host(context, ao_normalization);
    ao_normalization = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_ao_normalization",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_ao_factor_device(context, ao_normalization, ao_num, device_id);

  qmckl_free_host(context, ao_normalization);
  ao_normalization = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

  return QMCKL_SUCCESS;
}




qmckl_exit_code
qmckl_trexio_read_mo_X_device(qmckl_context_device context, trexio_t* const file, int device_id)
{
  assert (context != (qmckl_context_device) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;
  int64_t ao_num = 0L;

  // TODO
  rc = qmckl_get_ao_basis_ao_num_device(context, &ao_num);
  if (rc != QMCKL_SUCCESS)
    return rc;

int64_t mo_num = 0L;

rcio = trexio_read_mo_num_64(file, &mo_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_mo_num",
                         trexio_string_of_error(rcio));
}

assert (mo_num > 0);
// TODO
rc = qmckl_set_mo_basis_mo_num_device(context, mo_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ao_num * mo_num * sizeof(double);

  double* mo_coef = (double*) qmckl_malloc_host(context, mem_info);

  if (mo_coef == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_mo_X_device",
                           NULL);
  }

  assert (mo_coef != NULL);

  rcio = trexio_read_mo_coefficient_64(file, mo_coef);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_mo_coefficient",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_mo_basis_coefficient_device(context, mo_coef, device_id);

  qmckl_free_host(context, mo_coef);
  mo_coef = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

return QMCKL_SUCCESS;
}




qmckl_exit_code
qmckl_trexio_read_device(const qmckl_context_device context, const char* file_name, const int64_t size_max, int device_id)
{
  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_exit_code rc;
  char file_name_new[size_max+1];
  strncpy(file_name_new, file_name, size_max+1);
  file_name_new[size_max] = '\0';

  trexio_t* file = qmckl_trexio_open_X(file_name_new, &rc);
  if (file == NULL) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_trexio_read_device",
                           trexio_string_of_error(rc));
  }

  assert (file != NULL);

  rc = qmckl_trexio_read_electron_X(context, file);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_trexio_read_device",
                           "Error reading electron");
  }

  rc = qmckl_trexio_read_nucleus_X_device(context, file, device_id);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_trexio_read_device",
                           "Error reading nucleus");
  }

  rc = qmckl_trexio_read_ao_X_device(context, file, device_id);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_trexio_read_device",
                           "Error reading AOs");
  }

  rc = qmckl_trexio_read_mo_X_device(context, file, device_id);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_trexio_read",
                           "Error reading MOs");
  }

  trexio_close(file);
  file = NULL;
  return rc;
}
