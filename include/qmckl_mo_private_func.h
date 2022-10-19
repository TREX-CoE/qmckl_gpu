#ifndef QMCKL_MO_HPF
#define QMCKL_MO_HPF

/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */

qmckl_exit_code qmckl_init_mo_basis(qmckl_context context);

/* When the basis set is completely entered, other data structures are */
/* computed to accelerate the calculations. */

qmckl_exit_code qmckl_finalize_mo_basis(qmckl_context context);

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_value" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value(qmckl_context context);

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl(qmckl_context context);

#endif
