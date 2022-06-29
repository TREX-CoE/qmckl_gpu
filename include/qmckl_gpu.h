// This is the header file meant to be included by the users.
// It contains prototypes for all GPU functions.


// AO

#ifdef HAVE_OPENMP_OFFLOAD
qmckl_get_ao_basis_ao_vgl_omp_offload (qmckl_context context,
                                       double* const ao_vgl,
                                       const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace_omp_offload (qmckl_context context,
                                               double* const ao_vgl,
                                               const int64_t size_max);
#endif
