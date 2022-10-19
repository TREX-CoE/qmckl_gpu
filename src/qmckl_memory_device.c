#include "../include/qmckl_memory_device.h"


void* qmckl_malloc_host(qmckl_context_device context, const qmckl_memory_info_struct info) {

  assert (qmckl_context_check_device(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct_device* const ctx = (qmckl_context_struct_device*) context;

  /* Allocate memory and zero it */
  void * pointer = aligned_alloc(64, ((info.size+64) >> 6) << 6 );
  if (pointer == NULL) {
    return NULL;
  }
  memset(pointer, 0, info.size);

  qmckl_lock_device(context);
  {
    /* If qmckl_memory_struct is full, reallocate a larger one */
    if (ctx->memory.n_allocated == ctx->memory.array_size) {
      const size_t old_size = ctx->memory.array_size;
      qmckl_memory_info_struct * new_array = realloc(ctx->memory.element,
                                                          2L * old_size *
                                                          sizeof(qmckl_memory_info_struct));
      if (new_array == NULL) {
        qmckl_unlock_device(context);
        free(pointer);
        return NULL;
      }

      memset( &(new_array[old_size]), 0, old_size * sizeof(qmckl_memory_info_struct) );
      ctx->memory.element = new_array;
      ctx->memory.array_size = 2L * old_size;
    }

    /* Find first NULL entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].size > (size_t) 0) {
      pos += (size_t) 1;
    }
    assert (ctx->memory.element[pos].size == (size_t) 0);

    /* Copy info at the new location */
    memcpy(&(ctx->memory.element[pos]), &info, sizeof(qmckl_memory_info_struct));
    ctx->memory.element[pos].pointer = pointer;
    ctx->memory.n_allocated += (size_t) 1;
  }
  qmckl_unlock_device(context);

  return pointer;
}


qmckl_exit_code qmckl_free_host(qmckl_context_device context, void * const ptr) {

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith_device(context,
                            QMCKL_INVALID_CONTEXT,
                            "qmckl_free_hst",
                            NULL);
  }

  if (ptr == NULL) {
    return qmckl_failwith_device(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_free_host",
                          "NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  qmckl_lock_device(context);
  {
    /* Find pointer in array of saved pointers */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock_device(context);
      return qmckl_failwith_device(context,
                            QMCKL_INVALID_ARG_2,
                            "qmckl_free_host",
                            "Pointer not found in context");
    }

    free(ptr);

    memset( &(ctx->memory.element[pos]), 0, sizeof(qmckl_memory_info_struct) );
    ctx->memory.n_allocated -= (size_t) 1;
  }
  qmckl_unlock_device(context);

  return QMCKL_SUCCESS;
}



void* qmckl_malloc_device(qmckl_context_device context, const qmckl_memory_info_struct info, int device_id) {

  assert (qmckl_context_check_device(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct_device* const ctx = (qmckl_context_struct_device*) context;

  /* Allocate memory and zero it */
  void * pointer = omp_target_alloc(info.size, device_id);
  if (pointer == NULL) {
    return NULL;
  }

  // TODO
  // Memset to 0 of size info.size
  // memset(pointer, 0, info.size);

  qmckl_lock_device(context);
  {
    /* If qmckl_memory_struct is full, reallocate a larger one */
    if (ctx->memory.n_allocated == ctx->memory.array_size) {
      const size_t old_size = ctx->memory.array_size;
      qmckl_memory_info_struct * new_array = realloc(ctx->memory.element,
                                                          2L * old_size *
                                                          sizeof(qmckl_memory_info_struct));
      if (new_array == NULL) {
        qmckl_unlock_device(context);
        free(pointer);
        return NULL;
      }

      memset( &(new_array[old_size]), 0, old_size * sizeof(qmckl_memory_info_struct) );
      ctx->memory.element = new_array;
      ctx->memory.array_size = 2L * old_size;
    }

    /* Find first NULL entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].size > (size_t) 0) {
      pos += (size_t) 1;
    }
    assert (ctx->memory.element[pos].size == (size_t) 0);

    /* Copy info at the new location */
    ctx->memory.element[pos].size    = info.size;
    ctx->memory.element[pos].pointer = pointer;
    ctx->memory.n_allocated += (size_t) 1;
  }
  qmckl_unlock_device(context);

  return pointer;
}


qmckl_exit_code qmckl_free_device(qmckl_context_device context, void * const ptr, int device_id) {

  if (qmckl_context_check_device(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith_device(context,
                                   QMCKL_INVALID_CONTEXT,
                                   "qmckl_free_device",
                                   NULL);
  }

  if (ptr == NULL) {
    return qmckl_failwith_device(context,
                                 QMCKL_INVALID_ARG_2,
                                 "qmckl_free_device",
                                 "NULL pointer");
  }

  qmckl_context_struct_device* const ctx = (qmckl_context_struct_device*) context;

  qmckl_lock_device(context);
  {
    /* Find pointer in array of saved pointers */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock_device(context);
      return qmckl_failwith_device(context,
                                   QMCKL_FAILURE,
                                   "qmckl_free_device",
                                   "Pointer not found in context");
    }

    omp_target_free(ptr, device_id);

    memset( &(ctx->memory.element[pos]), 0, sizeof(qmckl_memory_info_struct) );
    ctx->memory.n_allocated -= (size_t) 1;
  }
  qmckl_unlock_device(context);

  return QMCKL_SUCCESS;
}
