#ifndef QMCKL_NUCLEUS_HPT
#define QMCKL_NUCLEUS_HPT
#include <stdbool.h>

/* Data structure */


typedef struct qmckl_nucleus_struct {
  int64_t      num;
  int64_t      repulsion_date;
  int64_t      nn_distance_date;
  int64_t      nn_distance_rescaled_date;
  int64_t      coord_date;
  qmckl_vector charge;
  qmckl_matrix coord;
  qmckl_matrix nn_distance;
  qmckl_matrix nn_distance_rescaled;
  double       repulsion;
  double       rescale_factor_kappa;
  int32_t      uninitialized;
  bool         provided;
} qmckl_nucleus_struct;

#endif
