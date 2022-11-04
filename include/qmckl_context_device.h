#pragma once

#include <qmckl.h>

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_blas_device.h"
#include "qmckl_device_types.h"
#include "qmckl_memory_device.h"


qmckl_context_device qmckl_context_create_device();

qmckl_exit_code
qmckl_context_destroy_device(const qmckl_context_device context, int device_id);
