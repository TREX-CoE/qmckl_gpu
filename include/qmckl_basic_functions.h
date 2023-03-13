#include <stdint.h>
#include <stdlib.h>
#include "qmckl_basic_type.h"

qmckl_context
qmckl_context_check (const qmckl_context context) ;

qmckl_exit_code
qmckl_context_touch (const qmckl_context context);


void qmckl_lock  (qmckl_context context);
void qmckl_unlock(qmckl_context context);

qmckl_context qmckl_context_create();

qmckl_exit_code
qmckl_context_destroy (const qmckl_context context);

const char*
qmckl_string_of_error (const qmckl_exit_code error);

qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message) ;

qmckl_exit_code qmckl_set_electron_num      (qmckl_context context, const int64_t up_num, const int64_t down_num);

qmckl_exit_code qmckl_set_electron_coord    (qmckl_context context, const char transp, const int64_t walk_num, const double* coord, const int64_t size_max);


qmckl_exit_code
qmckl_set_nucleus_num(qmckl_context context,
                      const int64_t num);

qmckl_exit_code
qmckl_set_nucleus_charge(qmckl_context context,
                         const double* charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_set_nucleus_coord(qmckl_context context,
                        const char transp,
                        const double* coord,
                        const int64_t size_max);

void* qmckl_malloc(qmckl_context context,
                   const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free(qmckl_context context,
                           void * const ptr);

qmckl_exit_code
qmckl_get_malloc_info(qmckl_context context,
                      const void* pointer,
                      qmckl_memory_info_struct* info);


typedef struct qmckl_memory_info_struct {
  size_t size;
  void*  pointer;
} qmckl_memory_info_struct;

static const qmckl_memory_info_struct qmckl_memory_info_struct_zero =
  {
   .size = (size_t) 0,
   .pointer = NULL
  };

typedef struct qmckl_memory_struct {
  size_t                    n_allocated;
  size_t                    array_size;
  qmckl_memory_info_struct* element;
} qmckl_memory_struct;




