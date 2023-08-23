// This file contains C functions that read arrays from n2.h, and that are
// meant to be called from the Fortran test

#include <inttypes.h>
#include "n2.h"

int64_t get_n2_nucl_num() { return n2_nucl_num; }

double* get_n2_charge() { return &(n2_charge[0]); }


int64_t get_n2_elec_up_num() { return n2_elec_up_num; }
int64_t get_n2_elec_dn_num() { return n2_elec_dn_num; }
int64_t get_n2_elec_num() { return n2_elec_num; }
int64_t get_n2_walk_num() { return n2_walk_num; }


double* get_n2_elec_coord() { return &(n2_elec_coord[0][0][0]); }
double* get_n2_nucl_coord() { return &(n2_nucl_coord[0][0]); }

int64_t get_n2_type_nucl_num() { return n2_type_nucl_num; }
int64_t get_n2_aord_num() { return n2_aord_num; }
int64_t get_n2_bord_num() { return n2_bord_num; }
int64_t get_n2_cord_num() { return n2_cord_num; }
int64_t get_n2_dim_cord_vec() { return n2_dim_cord_vec; }

int64_t* get_n2_type_nucl_vector() { return &(n2_type_nucl_vector[0]); }

double* get_n2_aord_vector() { return &(n2_aord_vector[0][0]); }
double* get_n2_bord_vector() { return &(n2_bord_vector[0]); }
double* get_n2_cord_vector() { return &(n2_cord_vector[0][0]); }
double* get_n2_cord_vector_full() { return &(n2_cord_vector_full[0][0]); }
double* get_n2_lkpm_of_cindex() { return &(n2_lkpm_of_cindex[0][0]); }

double* get_n2_rescale_factor_en() { return &(n2_rescale_factor_en[0]); }
double get_n2_rescale_factor_ee() { return n2_rescale_factor_ee; }
