#pragma once

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <pthread.h>

#include "qmckl_error_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_numprec_private_type.h"
#include "qmckl_point_private_type.h"
#include "qmckl_nucleus_private_type.h"
#include "qmckl_electron_private_type.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_mo_private_type.h"
#include "qmckl_jastrow_private_type.h"
#include "qmckl_determinant_private_type.h"
#include "qmckl_local_energy_private_type.h"
#include "qmckl_point_private_func.h"
#include "qmckl_nucleus_private_func.h"
#include "qmckl_electron_private_func.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_func.h"
#include "qmckl_determinant_private_func.h"
#include "qmckl_local_energy_private_func.h"

/* Data structure */


typedef struct qmckl_context_struct {
  /* -- State of the library -- */

  /* Validity checking */
  uint64_t                 tag;

  /* Numerical precision */
  qmckl_numprec_struct     numprec;

  /* Thread lock */
  int                      lock_count;
  pthread_mutex_t          mutex;

  /* Error handling */
  qmckl_error_struct       error;

  /* Memory allocation */
  qmckl_memory_struct      memory;

  /* Current date */
  uint64_t                 date;

  /* Points */
  qmckl_point_struct         point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct       nucleus;
  qmckl_electron_struct      electron;
  qmckl_ao_basis_struct      ao_basis;
  qmckl_mo_basis_struct      mo_basis;
  qmckl_jastrow_struct       jastrow;
  qmckl_determinant_struct   det;
  qmckl_local_energy_struct  local_energy;

  /* To be implemented:
  */

} qmckl_context_struct;



/* A tag is used internally to check if the memory domain pointed */
/* by a pointer is a valid context. This allows to check that even if */
/* the pointer associated with a context is non-null, we can still */
/* verify that it points to the expected data structure. */


#define VALID_TAG   0xBEEFFACE
#define INVALID_TAG 0xDEADBEEF
