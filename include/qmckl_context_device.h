#include <qmckl.h>

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_device_types.h"
#include "qmckl_memory_device.h"
#include "qmckl_blas_device.h"


qmckl_context_device qmckl_context_create_device();

qmckl_exit_code
qmckl_context_destroy_device (const qmckl_context_device context);

qmckl_context_device qmckl_context_check_device(const qmckl_context_device context);

qmckl_exit_code
qmckl_set_error_device(qmckl_context_device context,
                       const qmckl_exit_code exit_code,
                       const char* function_name,
                       const char* message);

qmckl_exit_code
qmckl_failwith_device(qmckl_context_device context,
                      const qmckl_exit_code exit_code,
                      const char* function,
                      const char* message);
