#pragma once

#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>

#include "qmckl_types.h"
#include "qmckl_types.h"

/* Error */
qmckl_exit_code_device
qmckl_failwith_device(qmckl_context_device context,
					  const qmckl_exit_code_device exit_code,
					  const char *function, const char *message);

/* Electron */
qmckl_exit_code_device
qmckl_set_electron_num_device(qmckl_context_device context,
							  const int64_t up_num, const int64_t down_num);
qmckl_exit_code_device qmckl_init_electron_device(qmckl_context_device context);
bool qmckl_nucleus_provided_device(qmckl_context_device context);
