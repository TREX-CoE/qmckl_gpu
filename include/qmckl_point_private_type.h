#ifndef QMCKL_POINT_HPT
#define QMCKL_POINT_HPT
#include <stdbool.h>
#include "qmckl_blas_private_type.h"

/* Data structure */


typedef struct qmckl_point_struct {
  int64_t      num;
  uint64_t     date;
  qmckl_matrix coord;
} qmckl_point_struct;

#endif
