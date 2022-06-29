#ifndef QMCKL_ELECTRON_HPT
#define QMCKL_ELECTRON_HPT
#include <stdbool.h>

/* Data structure */


typedef struct qmckl_electron_struct {
  int64_t        num;
  int64_t        up_num;
  int64_t        down_num;
  int64_t        walk_num;
  double         rescale_factor_kappa_ee;
  double         rescale_factor_kappa_en;
  int64_t        coord_new_date;
  int64_t        ee_distance_date;
  int64_t        en_distance_date;
  int64_t        ee_pot_date;
  int64_t        en_pot_date;
  int64_t        ee_distance_rescaled_date;
  int64_t        ee_distance_rescaled_deriv_e_date;
  int64_t        en_distance_rescaled_date;
  int64_t        en_distance_rescaled_deriv_e_date;
  qmckl_matrix   coord_new;
  qmckl_matrix   coord_old;
  double*        ee_distance;
  double*        en_distance;
  double*        ee_pot;
  double*        en_pot;
  double*        ee_distance_rescaled;
  double*        ee_distance_rescaled_deriv_e;
  double*        en_distance_rescaled;
  double*        en_distance_rescaled_deriv_e;
  int32_t        uninitialized;
  bool           provided;

  #ifdef HAVE_DEVICE_POINTERS
    double*        ee_distance_device;
    double*        en_distance_device;
    double*        ee_pot_device;
    double*        en_pot_device;
    double*        ee_distance_rescaled_device;
    double*        ee_distance_rescaled_deriv_e_device;
    double*        en_distance_rescaled_device;
    double*        en_distance_rescaled_deriv_e_device;
  #endif
} qmckl_electron_struct;

#endif
