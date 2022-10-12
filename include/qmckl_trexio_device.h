#include <stdint.h>
#include <trexio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include <stdio.h>

#include "qmckl.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"

#include "qmckl_device_types.h"

qmckl_exit_code
qmckl_trexio_read_nucleus_X_device(qmckl_context_device context, trexio_t* const file, int device_id);

qmckl_exit_code
qmckl_trexio_read_ao_X_device(qmckl_context context, trexio_t* const file, int device_id);

qmckl_exit_code
qmckl_trexio_read_mo_X_device(qmckl_context_device context, trexio_t* const file, int device_id);

qmckl_exit_code
qmckl_trexio_read_device(const qmckl_context_device context, const char* file_name, const int64_t size_max, int device_id);
