#pragma once

#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>

#include "qmckl_types.h"
#include "qmckl_context.h"

/* Error */
qmckl_exit_code_device
qmckl_failwith_device(qmckl_context_device context,
					  const qmckl_exit_code_device exit_code,
					  const char *function, const char *message);
