#ifndef QMCKL_ELECTRON_HPT
#define QMCKL_ELECTRON_HPT
#include "qmckl_point_private_type.h"
#include <stdbool.h>

/* Data structure */

typedef struct qmckl_walker_struct {
  int64_t num;
  qmckl_point_struct point;
} qmckl_walker;

typedef struct qmckl_electron_struct {
  int64_t num;
  int64_t up_num;
  int64_t down_num;
  qmckl_walker walker;
  qmckl_walker walker_old;
  uint64_t ee_distance_date;
  uint64_t en_distance_date;
  uint64_t ee_potential_date;
  uint64_t en_potential_date;
  double *ee_distance;
  double *en_distance;
  double *ee_potential;
  double *en_potential;
  int32_t uninitialized;
  bool provided;
} qmckl_electron_struct;

#endif
