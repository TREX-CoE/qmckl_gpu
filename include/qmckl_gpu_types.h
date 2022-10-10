// This file contains prototype for the GPU/device variants of the qmckl context structures.

#pragma once

#include <qmckl.h> // Needed of the datatypes that will be reused in GPU structures


//**********
// BLAS datatype (matrices, ...)
//**********



//**********
// Main datatype (qmckl_context_device)
//**********


typedef struct qmckl_context_struct_device {
  /* -- State of the library -- */

  /* Validity checking */
  uint64_t                 tag;

  /* Numerical precision */
  qmckl_numprec_struct              numprec;

  /* Thread lock */
  int                               lock_count;
  pthread_mutex_t                   mutex;

  /* Error handling */
  qmckl_error_struct                error;

  /* Memory allocation */
  qmckl_memory_struct_device        memory;

  /* Current date */
  uint64_t                          date;

  /* Points */
  qmckl_point_struct_devuce         point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct_device       nucleus;
  qmckl_electron_struct_device      electron;
  qmckl_ao_basis_struct_device      ao_basis;
  qmckl_mo_basis_struct             mo_basis; // TODO
  qmckl_jastrow_struct              jastrow; // TODO
  qmckl_determinant_struct          det;
  qmckl_local_energy_struct         local_energy;


} qmckl_context_struct_device;

typedef int64_t qmckl_context_device;
