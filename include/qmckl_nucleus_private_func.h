#ifndef QMCKL_NUCLEUS_HPF
#define QMCKL_NUCLEUS_HPF



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_nucleus(qmckl_context context);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_nn_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_nn_distance (
          const qmckl_context context,
          const int64_t nucl_num,
          const double* coord,
          double* const nn_distance );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_nn_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_compute_nn_distance_rescaled (
                                                    const qmckl_context context,
          const int64_t nucl_num,
          const double  rescale_factor_kappa,
          const double* coord,
          double* const nn_distance_rescaled );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_nucleus_repulsion(qmckl_context context);

qmckl_exit_code qmckl_compute_nucleus_repulsion (
     const qmckl_context context,
     const int64_t nucl_num,
     const double* charge,
     const double* nn_distance,
     double* energy
  );

#endif
