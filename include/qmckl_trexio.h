#pragma once

// This file contains functions prototypes for functions initializing the
// context from a TREXIO file

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <trexio.h>

#include "qmckl_types.h"
#include "qmckl_basic_functions.h"
#include "qmckl_nucleus.h"
#include "qmckl_electron.h"
#include "qmckl_context.h"
#include "qmckl_ao.h"
#include "qmckl_mo.h"

// Prototype for standard QMCkl function
trexio_t *qmckl_trexio_open_X_device(char *file_name,
									 qmckl_exit_code_device *rc);

qmckl_exit_code_device qmckl_trexio_read_device(qmckl_context_device context,
												char *file_name,
												int64_t size_max);

//**********
// CONTEXT FILL
//**********

qmckl_exit_code_device
qmckl_trexio_read_nucleus_X_device(qmckl_context_device context,
								   trexio_t *file);
qmckl_exit_code_device
qmckl_trexio_read_ao_X_device(qmckl_context_device context, trexio_t *file);
qmckl_exit_code_device
qmckl_trexio_read_mo_X_device(qmckl_context_device context, trexio_t *file);
qmckl_exit_code_device qmckl_trexio_read_device(qmckl_context_device context,
												char *file_name,
												int64_t size_max);
