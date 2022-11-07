#ifndef QMCKL_ELECTRON_HPT
#define QMCKL_ELECTRON_HPT
#include <stdbool.h>
#include "qmckl_point_private_type.h"

/* Data structure */


typedef struct qmckl_walker_struct {
  int64_t num;
  qmckl_point_struct point;
} qmckl_walker;
  
typedef struct qmckl_electron_struct {
  int64_t        num;
  int64_t        up_num;
  int64_t        down_num;
  qmckl_walker   walker;
  qmckl_walker   walker_old;
  double         rescale_factor_kappa_ee;
  double         rescale_factor_kappa_en;
  uint64_t        ee_distance_date;
  uint64_t        en_distance_date;
  uint64_t        ee_potential_date;
  uint64_t        en_potential_date;
  uint64_t        ee_distance_rescaled_date;
  uint64_t        ee_distance_rescaled_deriv_e_date;
  uint64_t        en_distance_rescaled_date;
  uint64_t        en_distance_rescaled_deriv_e_date;
  double*        ee_distance;
  double*        en_distance;
  double*        ee_potential;
  double*        en_potential;
  double*        ee_distance_rescaled;
  double*        ee_distance_rescaled_deriv_e;
  double*        en_distance_rescaled;
  double*        en_distance_rescaled_deriv_e;
  int32_t        uninitialized;
  bool           provided;

  #ifdef HAVE_DEVICE_POINTERS
  double*        ee_distance_device;
  double*        en_distance_device;
  double*        ee_potential_device;
  double*        en_potential_device;
  double*        ee_distance_rescaled_device;
  double*        ee_distance_rescaled_deriv_e_device;
  double*        en_distance_rescaled_device;
  double*        en_distance_rescaled_deriv_e_device;
  #endif
} qmckl_electron_struct;

#endif
