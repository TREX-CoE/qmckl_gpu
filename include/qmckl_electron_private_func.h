#ifndef QMCKL_ELECTRON_HPF
#define QMCKL_ELECTRON_HPF



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_electron(qmckl_context context);

/* CPU */


qmckl_exit_code qmckl_provide_ee_distance(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_ee_distance_device(qmckl_context context, int device_pointers);
#endif

qmckl_exit_code qmckl_compute_ee_distance (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance );

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_ee_distance_device (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance,
          int device_id );
#endif

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ee_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_kappa_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ee_distance_rescaled_deriv_e(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_distance_rescaled_deriv_e (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_kappa_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_deriv_e );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ee_potential(qmckl_context context);

/* CPU */


qmckl_exit_code qmckl_provide_en_distance(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_en_distance_device(qmckl_context context, int device_id);
#endif

qmckl_exit_code qmckl_compute_en_distance (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance );

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_en_distance_device (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance,
          int device_id
 );
#endif

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_en_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_compute_en_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const double  rescale_factor_kappa_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_en_distance_rescaled_deriv_e(qmckl_context context);

qmckl_exit_code qmckl_compute_en_distance_rescaled_deriv_e (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const double  rescale_factor_kappa_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_deriv_e );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_en_potential(qmckl_context context);

#endif
