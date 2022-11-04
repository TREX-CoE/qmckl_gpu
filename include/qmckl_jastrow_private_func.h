#ifndef QMCKL_JASTROW_HPF
#define QMCKL_JASTROW_HPF

/* Finalize basis (host & device) */

/*  When the required information is completely entered, other data structures are */
/*  computed to accelerate the calculations. The intermediates factors */
/*  are precontracted using BLAS LEVEL 3 operations for an optimal flop count. */


qmckl_exit_code qmckl_finalize_jastrow(qmckl_context context);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_finalize_jastrow_device(qmckl_context context);
#endif

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_asymp_jasb(qmckl_context context);



/* #   #+CALL: generate_c_header(table=qmckl_asymp_jasb_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_asymp_jasb (
      const qmckl_context context,
      const int64_t bord_num,
      const double* bord_vector,
      const double rescale_factor_kappa_ee,
      double* const asymp_jasb );

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_factor_ee(qmckl_context context);



/* #   #+CALL: generate_c_header(table=qmckl_factor_ee_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_ee (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t up_num,
      const int64_t bord_num,
      const double* bord_vector,
      const double* ee_distance_rescaled,
      const double* asymp_jasb,
      double* const factor_ee );

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_factor_ee_deriv_e(qmckl_context context);



/* #   #+CALL: generate_c_header(table=qmckl_factor_ee_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */



qmckl_exit_code qmckl_compute_factor_ee_deriv_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t up_num,
      const int64_t bord_num,
      const double* bord_vector,
      const double* ee_distance_rescaled,
      const double* ee_distance_rescaled_deriv_e,
      double* const factor_ee_deriv_e );

qmckl_exit_code qmckl_compute_factor_ee_deriv_e_hpc (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t up_num,
  const int64_t bord_num,
  const double* bord_vector,
  const double* ee_distance_rescaled,
  const double* ee_distance_rescaled_deriv_e,
  double* const factor_ee_deriv_e );

qmckl_exit_code qmckl_compute_factor_ee_deriv_e_doc (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t up_num,
  const int64_t bord_num,
  const double* bord_vector,
  const double* ee_distance_rescaled,
  const double* ee_distance_rescaled_deriv_e,
  double* const factor_ee_deriv_e );

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_factor_en(qmckl_context context);




/* #   #+CALL: generate_c_header(table=qmckl_factor_en_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_en (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const int64_t aord_num,
      const double* aord_vector,
      const double* en_distance_rescaled,
      double* const factor_en );

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_factor_en_deriv_e(qmckl_context context);



/* #   #+CALL: generate_c_header(table=qmckl_factor_en_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_en_deriv_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const int64_t aord_num,
      const double* aord_vector,
      const double* en_distance_rescaled,
      const double* en_distance_rescaled_deriv_e,
      double* const factor_en_deriv_e );

/* CPU */


qmckl_exit_code qmckl_provide_een_rescaled_e(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_e_device(qmckl_context context, int device_id);
#endif



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_rescaled_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_een_rescaled_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* ee_distance,
      double* const een_rescaled_e );



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_rescaled_e_args,rettyp=get_value("CRetType"),fname="qmckl_compute_een_rescaled_e_doc") */


qmckl_exit_code qmckl_compute_een_rescaled_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

qmckl_exit_code qmckl_compute_een_rescaled_e_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

qmckl_exit_code qmckl_compute_een_rescaled_e_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

/* CPU */


qmckl_exit_code qmckl_provide_een_rescaled_e_deriv_e(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_e_deriv_e_device(qmckl_context context, int device_id);
#endif



/* #  #+CALL: generate_c_header(table=qmckl_factor_een_rescaled_e_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een_rescaled_e_deriv_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* coord_ee,
      const double* ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_e_deriv_e );

/* Device pointers */



#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_factor_een_rescaled_e_deriv_e_device (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t cord_num,
         const double rescale_factor_kappa_ee,
         const double* coord_new,
         const double* ee_distance,
         const double* een_rescaled_e,
         double* const een_rescaled_e_deriv_e,
         int device_id
);
#endif

/* CPU */


qmckl_exit_code qmckl_provide_een_rescaled_n(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_n_device(qmckl_context context, int device_id);
#endif



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_rescaled_n_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_een_rescaled_n (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_en,
      const double* en_distance,
      double* const een_rescaled_n );

qmckl_exit_code qmckl_compute_een_rescaled_n_device (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_en,
      const double* en_distance,
      double* const een_rescaled_n,
      int device_id );

/* CPU */


qmckl_exit_code qmckl_provide_een_rescaled_n_deriv_e(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_n_deriv_e_device(qmckl_context context, int device_id);
#endif



/* #   #+CALL: generate_c_header(table=qmckl_compute_factor_een_rescaled_n_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een_rescaled_n_deriv_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_en,
      const double* coord_ee,
      const double* coord_en,
      const double* en_distance,
      const double* een_rescaled_n,
      double* const een_rescaled_n_deriv_e );

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
    qmckl_exit_code qmckl_compute_factor_een_rescaled_n_deriv_e_device (
          qmckl_context context,
          int64_t walk_num,
          int64_t elec_num,
          int64_t nucl_num,
          int64_t cord_num,
          double rescale_factor_kappa_en,
          double* coord_ee,
          double* coord_en,
          double* en_distance,
          double* een_rescaled_n,
          double* een_rescaled_n_deriv_e,
          int device_id );
#endif

/* CPU & offload version */


qmckl_exit_code qmckl_provide_dim_cord_vect(qmckl_context context);
qmckl_exit_code qmckl_provide_cord_vect_full(qmckl_context context);
qmckl_exit_code qmckl_provide_lkpm_combined_index(qmckl_context context);
qmckl_exit_code qmckl_provide_tmp_c(qmckl_context context);
qmckl_exit_code qmckl_provide_dtmp_c(qmckl_context context);

/* Device pointers version */


qmckl_exit_code qmckl_provide_cord_vect_full_device(qmckl_context context, int device_id);
qmckl_exit_code qmckl_provide_lkpm_combined_index_device(qmckl_context context, int device_id);
qmckl_exit_code qmckl_provide_tmp_c_device(qmckl_context context, int device_id);
qmckl_exit_code qmckl_provide_dtmp_c_device(qmckl_context context, int device_id);



/* #   #+CALL: generate_c_header(table=qmckl_factor_dim_cord_vect_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_dim_cord_vect (
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_cord_vect );




/* #   #+CALL: generate_c_header(table=qmckl_factor_cord_vect_full_args,rettyp=get_value("CRetType"),fname="qmckl_compute_cord_vect_full_doc") */


qmckl_exit_code qmckl_compute_cord_vect_full (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_cord_vect,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* cord_vector,
      double* const cord_vect_full );

qmckl_exit_code qmckl_compute_cord_vect_full_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_cord_vect,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* cord_vector,
      double* const cord_vect_full );

qmckl_exit_code qmckl_compute_cord_vect_full_hpc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_cord_vect,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* cord_vector,
      double* const cord_vect_full );

/* Device pointers version */


qmckl_exit_code qmckl_compute_cord_vect_full_device (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_cord_vect,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* cord_vector,
	  double* const cord_vect_full,
	  int device_id
);



/* #   #+CALL: generate_c_header(table=qmckl_factor_lkpm_combined_index_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_cord_vect,
      int64_t* const lkpm_combined_index );

/* Device pointers version */


qmckl_exit_code qmckl_compute_lkpm_combined_index_device (
        const qmckl_context context,
        const int64_t cord_num,
        const int64_t dim_cord_vect,
        int64_t* const lkpm_combined_index,
        int device_id
);



/* #   #+CALL: generate_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c") */


qmckl_exit_code qmckl_compute_tmp_c (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e,
      const double* een_rescaled_n,
      double* const tmp_c );

qmckl_exit_code qmckl_compute_tmp_c_doc (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );



/* #   #+CALL: generate_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c_doc") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_compute_tmp_c_doc (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );



/* #   #+CALL: generate_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c_hpc") */

/*     #+RESULTS: */


qmckl_exit_code qmckl_compute_tmp_c_hpc (const qmckl_context context,
                                         const int64_t cord_num,
                                         const int64_t elec_num,
                                         const int64_t nucl_num,
                                         const int64_t walk_num,
                                         const double* een_rescaled_e,
                                         const double* een_rescaled_n,
                                         double* const tmp_c );

#ifdef HAVE_OPENACC_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_acc_offload (const qmckl_context context,
                                 const int64_t cord_num,
                                 const int64_t elec_num,
                                 const int64_t nucl_num,
                                 const int64_t walk_num,
                                 const double* een_rescaled_e,
                                 const double* een_rescaled_n,
                                 double* const tmp_c );
#endif

#ifdef HAVE_OPENMP_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_omp_offload (const qmckl_context context,
                                 const int64_t cord_num,
                                 const int64_t elec_num,
                                 const int64_t nucl_num,
                                 const int64_t walk_num,
                                 const double* een_rescaled_e,
                                 const double* een_rescaled_n,
                                 double* const tmp_c );
#endif

#ifdef HAVE_CUBLAS_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_cublas_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e,
      const double* een_rescaled_n,
      double* const tmp_c );
#endif

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_compute_tmp_c_device (const qmckl_context context,
                            const int64_t cord_num,
                            const int64_t elec_num,
                            const int64_t nucl_num,
                            const int64_t walk_num,
                            const double* een_rescaled_e,
                            const double* een_rescaled_n,
                            double* const tmp_c,

                            int device_id );
#endif

/* Doc */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_dtmp_c */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*     #+NAME: qmckl_factor_dtmp_c_args */
/*     | Variable                 | Type                                                             | In/Out | Description                                   | */
/*     |--------------------------+------------------------------------------------------------------+--------+-----------------------------------------------| */
/*     | ~context~                | ~qmckl_context~                                                  | in     | Global state                                  | */
/*     | ~cord_num~               | ~int64_t~                                                        | in     | Order of polynomials                          | */
/*     | ~elec_num~               | ~int64_t~                                                        | in     | Number of electrons                           | */
/*     | ~nucl_num~               | ~int64_t~                                                        | in     | Number of nucleii                             | */
/*     | ~walk_num~               | ~int64_t~                                                        | in     | Number of walkers                             | */
/*     | ~een_rescaled_e_deriv_e~ | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~            | in     | Electron-electron rescaled factor derivatives | */
/*     | ~een_rescaled_n~         | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled factor              | */
/*     | ~dtmp_c~                 | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | out    | vector of non-zero coefficients               | */



qmckl_exit_code
qmckl_compute_dtmp_c (const qmckl_context context,
                      const int64_t cord_num,
                      const int64_t elec_num,
                      const int64_t nucl_num,
                      const int64_t walk_num,
                      const double* een_rescaled_e_deriv_e,
                      const double* een_rescaled_n,
                      double* const dtmp_c );

qmckl_exit_code qmckl_compute_dtmp_c_doc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );

qmckl_exit_code qmckl_compute_dtmp_c_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );

#ifdef HAVE_OPENACC_OFFLOAD
    qmckl_exit_code qmckl_compute_dtmp_c_acc_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );
#endif

#ifdef HAVE_OPENMP_OFFLOAD
    qmckl_exit_code qmckl_compute_dtmp_c_omp_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );
#endif



/* #+RESULTS: */

#ifdef HAVE_CUBLAS_OFFLOAD
    qmckl_exit_code qmckl_compute_dtmp_c_cublas_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );
#endif

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
    qmckl_exit_code qmckl_compute_dtmp_c_omp_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c );
#endif

/* Provide                                                        :noexport: */

qmckl_exit_code qmckl_provide_factor_een(qmckl_context context);



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een_naive (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_cord_vect,
  const double* cord_vect_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_e,
  const double* een_rescaled_n,
  double* const factor_een );



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_cord_vect,
  const double* cord_vect_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_e,
  const double* een_rescaled_n,
  double* const factor_een );

/* CPU */


qmckl_exit_code qmckl_provide_factor_een_deriv_e(qmckl_context context);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_factor_een_deriv_e_device(qmckl_context context, int device_id);
#endif



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_deriv_e_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een_deriv_e_naive (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_cord_vect,
  const double* cord_vect_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_e,
  const double* een_rescaled_n,
  const double* een_rescaled_e_deriv_e,
  const double* een_rescaled_n_deriv_e,
  double* const factor_een_deriv_e );



/* #   #+CALL: generate_c_header(table=qmckl_factor_een_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_factor_een_deriv_e (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_cord_vect,
  const double* cord_vect_full,
  const int64_t* lkpm_combined_index,
  const double* tmp_c,
  const double* dtmp_c,
  const double* een_rescaled_n,
  const double* een_rescaled_n_deriv_e,
  double* const factor_een_deriv_e );

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
    qmckl_exit_code qmckl_compute_factor_een_deriv_e_device (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_cord_vect,
      const double* cord_vect_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_deriv_e,
      double* const factor_een_deriv_e,
	  int device_id );
#endif

#endif
