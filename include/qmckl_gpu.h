// This is the header file meant to be included by the users.
// It contains prototypes for all GPU functions, and definition of
// _device variants fo the structure datatypes.

#pragma once

#include <qmckl.h>

//**********
// GPU/DEVICE DATATYPES
//**********

#include "qmckl_gpu_datatypes.h"

// BLAS types (matrices, ...)

typedef struct qmckl_vector_device {
  double* restrict data;
  double* restrict data_device;

  int64_t size;
} qmckl_vector_device;

typedef struct qmckl_matrix_device {
  double* restrict data;
  double* restrict data_device;

  int64_t size[2];
} qmckl_matrix_device;


// Context subtypes
typedef struct qmckl_context_struct_device {
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
  qmckl_point_struct_device         point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct_device       nucleus;
  qmckl_electron_struct_device      electron;
  qmckl_ao_basis_struct_device      ao_basis;
  qmckl_mo_basis_struct             mo_basis; // TODO
  qmckl_jastrow_struct              jastrow; // TODO
  qmckl_determinant_struct          det;
  qmckl_local_energy_struct         local_energy;

} qmckl_context_struct;

// Main context type


//**********
// FUNCTION HEADERS
//**********

// AO
#include "qmckl_ao_gpu.h"
