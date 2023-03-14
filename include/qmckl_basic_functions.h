#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>
#include "qmckl_basic_type.h"

/* Context */
qmckl_context
qmckl_context_check (const qmckl_context context) ;

qmckl_exit_code
qmckl_context_touch (const qmckl_context context);


void qmckl_lock  (qmckl_context context);
void qmckl_unlock(qmckl_context context);

qmckl_context qmckl_context_create();

qmckl_exit_code
qmckl_context_destroy (const qmckl_context context);

/* Error */
const char*
qmckl_string_of_error (const qmckl_exit_code error);

qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message) ;

/* Electron */
qmckl_exit_code qmckl_set_electron_num      (qmckl_context context, const int64_t up_num, const int64_t down_num);

qmckl_exit_code qmckl_set_electron_coord    (qmckl_context context, const char transp, const int64_t walk_num, const double* coord, const int64_t size_max);


qmckl_exit_code qmckl_init_electron(qmckl_context context);

/* Nucleus */

qmckl_exit_code qmckl_init_nucleus(qmckl_context context);

qmckl_exit_code
qmckl_get_nucleus_num (const qmckl_context context, int64_t* const numi);

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

/* Memory */


qmckl_exit_code
qmckl_set_error(qmckl_context context,
                const qmckl_exit_code exit_code,
                const char* function_name,
                const char* message);

void* qmckl_malloc(qmckl_context context,
                   const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free(qmckl_context context,
                           void * const ptr);

qmckl_exit_code
qmckl_get_malloc_info(qmckl_context context,
                      const void* pointer,
                      qmckl_memory_info_struct* info);



/* Point */

qmckl_exit_code
qmckl_set_point (qmckl_context context,
                 const char transp,
                 const int64_t num,
                 const double* coord,
                 const int64_t size_max);

qmckl_exit_code qmckl_init_point(qmckl_context context);


/* AO */

qmckl_exit_code qmckl_init_ao_basis(qmckl_context context);
qmckl_exit_code
qmckl_get_ao_basis_ao_num (const qmckl_context context,
                           int64_t* const ao_num);
qmckl_exit_code
qmckl_get_ao_basis_shell_num (const qmckl_context context,
                              int64_t* const shell_num);
qmckl_exit_code
qmckl_get_ao_basis_prim_num (const qmckl_context context,
                             int64_t* const prim_num);
qmckl_exit_code
qmckl_get_ao_basis_type (const qmckl_context context,
                         char* const basis_type);

/* MO */

qmckl_exit_code qmckl_init_mo_basis(qmckl_context context);
qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num);

/* Determinant */

qmckl_exit_code qmckl_init_determinant(qmckl_context context);


/* Jastrow */

qmckl_exit_code qmckl_init_jastrow(qmckl_context context);


/* BLAS */

qmckl_vector
qmckl_vector_alloc( qmckl_context context,
                    const int64_t size);


qmckl_exit_code
qmckl_vector_free( qmckl_context context,
                   qmckl_vector* vector);


qmckl_matrix
qmckl_matrix_alloc( qmckl_context context,
                    const int64_t size1,
                    const int64_t size2);


qmckl_exit_code
qmckl_matrix_free( qmckl_context context,
                   qmckl_matrix* matrix);


qmckl_exit_code
qmckl_vector_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_vector* vector_out);



qmckl_exit_code
qmckl_matrix_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_matrix* matrix);


qmckl_exit_code
qmckl_transpose (qmckl_context context,
                 const qmckl_matrix A,
                 qmckl_matrix At );


qmckl_vector
qmckl_vector_of_matrix(const qmckl_matrix matrix);


qmckl_matrix
qmckl_matrix_of_vector(const qmckl_vector vector,
                       const int64_t size1,
                       const int64_t size2);

