#ifndef QMCKL_MO_HPT
#define QMCKL_MO_HPT

#include <stdbool.h>

/* Data structure */


typedef struct qmckl_mo_basis_struct {
  int64_t  mo_num;
  double * restrict coefficient;
  double * restrict coefficient_t;

  double * restrict mo_vgl;
  double * restrict mo_value;
  uint64_t  mo_vgl_date;
  uint64_t  mo_value_date;

  int32_t   uninitialized;
  bool      provided;
} qmckl_mo_basis_struct;

#endif
