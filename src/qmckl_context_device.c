// This file will contains general functions accepting a qmckl_contect_device as
// argument, hence why they need a device alternative

#include "../include/qmckl_context_device.h"


//**********
// CONTEXT CREATE/DESTROY
//**********

qmckl_exit_code qmckl_init_point_device(qmckl_context_device context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL);

  memset(&(ctx->point), 0, sizeof(qmckl_point_struct));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_init_electron_device(qmckl_context_device context) {

  if (qmckl_context_check((qmckl_context) context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL);

  ctx->electron.uninitialized = (1 << 1) - 1;

  /* Default values */
  ctx->electron.rescale_factor_kappa_ee = 1.0;
  ctx->electron.rescale_factor_kappa_en = 1.0;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_init_nucleus_device(qmckl_context_device context) {

  if (qmckl_context_check((qmckl_context) context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL);

  ctx->nucleus.uninitialized = (1 << 3) - 1;

  /* Default values */
  ctx->nucleus.rescale_factor_kappa = 1.0;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_init_ao_basis_device(qmckl_context_device context) {

  if (qmckl_context_check((qmckl_context) context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context, QMCKL_INVALID_CONTEXT,
                          "qmckl_init_ao_basis_device", NULL);
  }

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL);

  ctx->ao_basis.uninitialized = (1 << 14) - 1;

  /* Default values */
  ctx->ao_basis.ao_cartesian = true;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_init_mo_basis_device(qmckl_context_device context) {

  if (qmckl_context_check((qmckl_context) context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL);

  ctx->mo_basis.uninitialized = (1 << 2) - 1;

  return QMCKL_SUCCESS;
}

qmckl_context_device qmckl_context_create_device() {

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)malloc(
          sizeof(qmckl_context_struct));

  if (ctx == NULL) {
    return QMCKL_NULL_CONTEXT;
  }

  /* Set all pointers and values to NULL */
  { memset(ctx, 0, sizeof(qmckl_context_struct)); }

  /* Initialize lock */
  {
    pthread_mutexattr_t attr;
    int rc;

    rc = pthread_mutexattr_init(&attr);
    assert(rc == 0);

#ifdef PTHREAD_MUTEX_RECURSIVE
    (void)pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
#endif

    rc = pthread_mutex_init(&(ctx->mutex), &attr);
    assert(rc == 0);

    (void)pthread_mutexattr_destroy(&attr);
  }

  /* Initialize data */
  {
    ctx->tag = VALID_TAG;

    const qmckl_context context = (qmckl_context)ctx;
    assert(qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

    qmckl_exit_code rc;

    ctx->numprec.precision = QMCKL_DEFAULT_PRECISION;
    ctx->numprec.range = QMCKL_DEFAULT_RANGE;

    // NOTE Keeping device variants for now, but everything could
    // probably be relplaced (including this function)
    rc = qmckl_init_point_device(context);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_init_electron_device(context);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_init_nucleus_device(context);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_init_ao_basis_device(context);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_init_mo_basis_device(context);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_init_determinant(context);
    assert(rc == QMCKL_SUCCESS);
  }

  /* Allocate qmckl_memory_struct */
  {
    const size_t size = 128L;
    qmckl_memory_info_struct *new_array =
        calloc(size, sizeof(qmckl_memory_info_struct));
    if (new_array == NULL) {
      return QMCKL_NULL_CONTEXT;
    }
    memset(&(new_array[0]), 0, size * sizeof(qmckl_memory_info_struct));

    ctx->memory.element = new_array;
    ctx->memory.array_size = size;
    ctx->memory.n_allocated = (size_t)0;
  }

  return (qmckl_context_device)ctx;
}

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context, int device_id) {

  const qmckl_context_device checked_context =
      qmckl_context_check((qmckl_context) context);
  if (checked_context == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  qmckl_context_struct *const ctx =
      (qmckl_context_struct *)context;
  assert(ctx != NULL); /* Shouldn't be possible because the context is valid */

  qmckl_lock((qmckl_context) context);
  {
    /* Memory: Remove all allocated data */
    // TODO This should free memory on device
    for (size_t pos = (size_t)0; pos < ctx->memory.array_size; ++pos) {
      if (ctx->memory.element[pos].pointer != NULL) {
        omp_target_free(ctx->memory.element[pos].pointer, device_id);
        memset(&(ctx->memory.element[pos]), 0,
               sizeof(qmckl_memory_info_struct));
        ctx->memory.n_allocated -= 1;
      }
    }
    assert(ctx->memory.n_allocated == (size_t)0);
    free(ctx->memory.element);
    ctx->memory.element = NULL;
    ctx->memory.array_size = (size_t)0;
  }
  qmckl_unlock((qmckl_context) context);

  ctx->tag = INVALID_TAG;

  const int rc_destroy = pthread_mutex_destroy(&(ctx->mutex));
  if (rc_destroy != 0) {
    /* DEBUG */
    fprintf(stderr, "qmckl_context_destroy: %s (count = %d)\n",
            strerror(rc_destroy), ctx->lock_count);
    abort();
  }

  free(ctx);

  return QMCKL_SUCCESS;
}
