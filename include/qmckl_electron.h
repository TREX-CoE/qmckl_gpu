#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_context.h"
#include "qmckl_memory.h"
#include "qmckl_blas.h"
#include "qmckl_point.h"
#include "qmckl_distance.h"

bool qmckl_electron_provided_device(qmckl_context_device context);

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
qmckl_exit_code_device
qmckl_set_electron_coord_device_from_host(qmckl_context_device context, char transp,
								int64_t walk_num, double *coord,
								int64_t size_max);


qmckl_exit_code_device
qmckl_get_electron_coord_device(const qmckl_context_device context,
								const char transp, double *const coord,
								const int64_t size_max);

qmckl_exit_code_device qmckl_compute_en_distance_device(
	const qmckl_context_device context, const int64_t point_num,
	const int64_t nucl_num, const double *elec_coord, const double *nucl_coord,
	double *const en_distance);

qmckl_exit_code_device qmckl_compute_ee_distance_device(
	const qmckl_context_device context, const int64_t elec_num,
	const int64_t walk_num, const double *coord, double *const ee_distance);
