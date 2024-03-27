// This file contains C functions that read arrays from chbrclf.h, and that are
// meant to be called from the Fortran test

#include <inttypes.h>
#include "chbrclf.h"

int64_t get_chbrclf_elec_up_num() { return chbrclf_elec_up_num; }

int64_t get_chbrclf_elec_dn_num() { return chbrclf_elec_dn_num; }

int64_t get_chbrclf_nucl_num() { return chbrclf_nucl_num; }

int64_t get_chbrclf_walk_num() { return chbrclf_walk_num; }

#include <stdio.h>
int64_t get_chbrclf_elec_num() { return chbrclf_elec_num; }

int64_t get_chbrclf_shell_num() { return chbrclf_shell_num; }

int64_t get_chbrclf_prim_num() { return chbrclf_prim_num; }

int64_t get_chbrclf_ao_num() { return chbrclf_ao_num; }

int64_t get_chbrclf_mo_num() { return chbrclf_mo_num; }

double *get_chbrclf_elec_coord() { return &(chbrclf_elec_coord[0][0][0]); }

double *get_chbrclf_nucl_charge() { return chbrclf_charge; }

double *get_chbrclf_nucl_coord() { return &(chbrclf_nucl_coord[0][0]); }

int64_t *get_chbrclf_nucleus_index() {
	return &(chbrclf_basis_nucleus_index[0]);
}

int64_t *get_chbrclf_nucleus_shell_num() {
	return &(chbrclf_basis_nucleus_shell_num[0]);
}

int32_t *get_chbrclf_shell_ang_mom() {
	return &(chbrclf_basis_shell_ang_mom[0]);
}

int64_t *get_chbrclf_shell_prim_num() {
	return &(chbrclf_basis_shell_prim_num[0]);
}

int64_t *get_chbrclf_shell_prim_index() {
	return &(chbrclf_basis_shell_prim_index[0]);
}

double *get_chbrclf_shell_factor() { return &(chbrclf_basis_shell_factor[0]); }

double *get_chbrclf_exponent() { return &(chbrclf_basis_exponent[0]); }

double *get_chbrclf_coefficient() { return &(chbrclf_basis_coefficient[0]); }

double *get_chbrclf_prim_factor() { return &(chbrclf_basis_prim_factor[0]); }

double *get_chbrclf_ao_factor() { return &(chbrclf_basis_ao_factor[0]); }

double *get_chbrclf_mo_coef() { return &(chbrclf_mo_coef[0]); }
