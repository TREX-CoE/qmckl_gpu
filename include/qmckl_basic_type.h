#include <stdbool.h>

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
  uint64_t       ee_distance_date;
  uint64_t       en_distance_date;
  uint64_t       ee_potential_date;
  uint64_t       en_potential_date;
  double*        ee_distance;
  double*        en_distance;
  double*        ee_potential;
  double*        en_potential;
  int32_t        uninitialized;
  bool           provided;
} qmckl_electron_struct;

typedef struct qmckl_context_struct {
  /* -- State of the library -- */

  /* Validity checking */
  uint64_t                 tag;

  /* Numerical precision */
  qmckl_numprec_struct     numprec;

  /* Thread lock */
  int                      lock_count;
  pthread_mutex_t          mutex;

  /* Error handling */
  qmckl_error_struct       error;

  /* Memory allocation */
  qmckl_memory_struct      memory;

  /* Current date */
  uint64_t                 date;

  /* Points */
  qmckl_point_struct         point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct       nucleus;
  qmckl_electron_struct      electron;
  qmckl_ao_basis_struct      ao_basis;
  qmckl_mo_basis_struct      mo_basis;
  qmckl_jastrow_struct       jastrow;
  qmckl_determinant_struct   det;
  qmckl_local_energy_struct  local_energy;

  /* To be implemented:
  */

  /* Pointer to implementation-specific data */

  void* qmckl_extra;

} qmckl_context_struct;

typedef struct qmckl_point_struct {
  int64_t      num;
  uint64_t     date;
  qmckl_matrix coord;
} qmckl_point_struct;


