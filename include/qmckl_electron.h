#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"
#include "qmckl_point.h"

qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context, int64_t up_num,
							  int64_t down_num);

qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context, int64_t up_num,
							  int64_t down_num);

qmckl_exit_code_device
qmckl_set_electron_coord_device(qmckl_context_device context, char transp,
								int64_t walk_num, double *coord,
								int64_t size_max);
