#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"

qmckl_exit_code_device qmckl_set_point_device(qmckl_context_device context,
											  char transp, int64_t num,
											  double *coord, int64_t size_max);

qmckl_exit_code_device
qmckl_get_point_device(const qmckl_context_device context, const char transp,
					   double *const coord, const int64_t size_max);
