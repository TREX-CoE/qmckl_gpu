#ifndef QMCKL_JASTROW_HPT
#define QMCKL_JASTROW_HPT
#include <stdbool.h>

/* Data structure */

typedef struct qmckl_jastrow_struct {
  int32_t uninitialized;
  int64_t aord_num;
  int64_t bord_num;
  int64_t cord_num;
  int64_t type_nucl_num;
  uint64_t asymp_jasb_date;
  uint64_t tmp_c_date;
  uint64_t dtmp_c_date;
  uint64_t factor_ee_date;
  uint64_t factor_en_date;
  uint64_t factor_een_date;
  uint64_t factor_ee_deriv_e_date;
  uint64_t factor_en_deriv_e_date;
  uint64_t factor_een_deriv_e_date;
  double rescale_factor_ee;
  double *rescale_factor_en;
  int64_t *type_nucl_vector;
  double *aord_vector;
  double *bord_vector;
  double *cord_vector;
  double *asymp_jasb;
  double *factor_ee;
  double *factor_en;
  double *factor_een;
  double *factor_ee_deriv_e;
  double *factor_en_deriv_e;
  double *factor_een_deriv_e;
  int64_t dim_cord_vect;
  uint64_t dim_cord_vect_date;
  double *cord_vect_full;
  uint64_t cord_vect_full_date;
  int64_t *lkpm_combined_index;
  uint64_t lkpm_combined_index_date;
  double *tmp_c;
  double *dtmp_c;
  uint64_t ee_distance_rescaled_date;
  uint64_t ee_distance_rescaled_deriv_e_date;
  uint64_t en_distance_rescaled_date;
  uint64_t en_distance_rescaled_deriv_e_date;
  double *ee_distance_rescaled;
  double *ee_distance_rescaled_deriv_e;
  double *en_distance_rescaled;
  double *en_distance_rescaled_deriv_e;
  double *een_rescaled_e;
  double *een_rescaled_n;
  uint64_t een_rescaled_e_date;
  uint64_t een_rescaled_n_date;
  double *een_rescaled_e_deriv_e;
  double *een_rescaled_n_deriv_e;
  uint64_t een_rescaled_e_deriv_e_date;
  uint64_t een_rescaled_n_deriv_e_date;
  bool provided;
  char *type;

#ifdef HAVE_HPC
  bool gpu_offload;
#endif
} qmckl_jastrow_struct;

#endif
