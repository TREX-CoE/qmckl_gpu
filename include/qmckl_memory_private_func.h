#ifndef QMCKL_MEMORY_HPF
#define QMCKL_MEMORY_HPF

void* qmckl_malloc(qmckl_context context,
                   const qmckl_memory_info_struct info);

#ifdef HAVE_DEVICE_POINTERS
void* qmckl_malloc_device(qmckl_context context,
                          const qmckl_memory_info_struct info,
	                        int device_id);
#endif

qmckl_exit_code qmckl_free(qmckl_context context,
                           void * const ptr);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_free_device(qmckl_context context,
                                  void * const ptr,
	                                int device_id);
#endif

/* End of files                                                     :noexport: */


#endif
