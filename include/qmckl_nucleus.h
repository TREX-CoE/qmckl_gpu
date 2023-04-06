
#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"

qmckl_exit_code_device
qmckl_get_nucleus_num_device(qmckl_context_device context, int64_t *num);
qmckl_exit_code_device
qmckl_get_nucleus_num_device(qmckl_context_device context, int64_t *num);

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num);

qmckl_exit_code_device
qmckl_set_nucleus_num_device(qmckl_context_device context, int64_t num);

qmckl_exit_code_device
qmckl_set_nucleus_charge_device(qmckl_context_device context, double *charge,
								int64_t size_max);
qmckl_exit_code_device
qmckl_set_nucleus_coord_device(qmckl_context_device context, char transp,
							   double *coord, int64_t size_max);

qmckl_exit_code_device
qmckl_finalize_nucleus_basis_hpc_device(qmckl_context_device context);
qmckl_exit_code_device
qmckl_finalize_nucleus_basis_device(qmckl_context_device context);
