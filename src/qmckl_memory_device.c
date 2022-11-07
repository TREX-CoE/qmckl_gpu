#include "../include/qmckl_memory_device.h"

//**********
// HOST MEMORY
//**********

void *qmckl_malloc_host(qmckl_context_device context,
                        const qmckl_memory_info_struct info) {
  qmckl_context context_base = (qmckl_context)context;
  return qmckl_malloc(context_base, info);
}

qmckl_exit_code qmckl_free_host(qmckl_context_device context, void *const ptr) {
  qmckl_context context_base = (qmckl_context)context;
  return qmckl_free(context_base, ptr);
}
