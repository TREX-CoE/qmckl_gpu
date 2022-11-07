#ifndef QMCKL_MEMORY_HPF
#define QMCKL_MEMORY_HPF

void *qmckl_malloc(qmckl_context context, const qmckl_memory_info_struct info);

#ifdef HAVE_DEVICE_POINTERS
void *qmckl_malloc_device(qmckl_context context,
                          const qmckl_memory_info_struct info, int device_id);
#endif

qmckl_exit_code qmckl_free(qmckl_context context, void *const ptr);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_free_device(qmckl_context context, void *const ptr,
                                  int device_id);
#endif

qmckl_exit_code qmckl_get_malloc_info(qmckl_context context,
                                      const void *pointer,
                                      qmckl_memory_info_struct *info);

qmckl_exit_code qmckl_get_malloc_info_device(qmckl_context context,
                                             const void *pointer,
                                             qmckl_memory_info_struct *info,
                                             const int64_t device_id);

/* End of files                                                     :noexport:
 */

#endif
