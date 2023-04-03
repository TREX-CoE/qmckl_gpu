#include "include/qmckl_basic_functions.h"

#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <errno.h>
#include <math.h>

/* Error */

qmckl_exit_code_device
qmckl_set_error_device(qmckl_context_device context,
					   const qmckl_exit_code_device exit_code,
					   const char *function_name, const char *message) {
	/* Passing a function name and a message is mandatory. */
	assert(function_name != NULL);
	assert(message != NULL);

	/* Exit codes are assumed valid. */
	assert(exit_code >= 0);
	assert(exit_code != QMCKL_SUCCESS_DEVICE);
	assert(exit_code < QMCKL_INVALID_EXIT_CODE_DEVICE);

	/* The context is assumed to exist. */
	assert(qmckl_context_check(context) != QMCKL_NULL_CONTEXT_DEVICE);

	qmckl_lock_device(context);
	{
		qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
		assert(ctx != NULL); /* Impossible because the context is valid. */

		ctx->error.exit_code = exit_code;
		strncpy(ctx->error.function, function_name, QMCKL_MAX_FUN_LEN_DEVICE - 1);
		strncpy(ctx->error.message, message, QMCKL_MAX_MSG_LEN_DEVICE - 1);
	}
	qmckl_unlock_device(context);

	return QMCKL_SUCCESS_DEVICE;
}

const char *qmckl_string_of_error_device(const qmckl_exit_code_device error) {
	switch (error) {}
	return "Unknown error";
}

qmckl_exit_code_device
qmckl_failwith_device(qmckl_context_device context,
					  const qmckl_exit_code_device exit_code,
					  const char *function, const char *message) {
	assert(exit_code > 0);
	assert(exit_code < QMCKL_INVALID_EXIT_CODE_DEVICE);
	assert(function != NULL);
	assert(strlen(function) < QMCKL_MAX_FUN_LEN_DEVICE);
	if (message != NULL) {
		assert(strlen(message) < QMCKL_MAX_MSG_LEN_DEVICE);
	}

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT_DEVICE)
		return QMCKL_INVALID_CONTEXT_DEVICE;

	if (message == NULL) {
		qmckl_exit_code_device rc = qmckl_set_error(
			context, exit_code, function, qmckl_string_of_error(exit_code));
		assert(rc == QMCKL_SUCCESS_DEVICE);
	} else {
		qmckl_exit_code_device rc =
			qmckl_set_error(context, exit_code, function, message);
		assert(rc == QMCKL_SUCCESS_DEVICE);
	}

	return exit_code;
}

/* Electron */


bool qmckl_electron_provided_device(const qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->electron.provided;
}

qmckl_exit_code_device
qmckl_get_electron_num_device(const qmckl_context_device context,
							  int64_t *const num) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return QMCKL_INVALID_CONTEXT_DEVICE;
	}

	if (num == NULL) {
		return qmckl_failwith(context, QMCKL_INVALID_ARG_2_DEVICE,
							  "qmckl_get_electron_num",
							  "num is a null pointer");
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	int32_t mask = 1 << 0;

	if ((ctx->electron.uninitialized & mask) != 0) {
		return QMCKL_NOT_PROVIDED;
	}

	assert(ctx->electron.num > (int64_t)0);
	*num = ctx->electron.num;
	return QMCKL_SUCCESS_DEVICE;
}

/* nucleus */


bool qmckl_nucleus_provided_device(qmckl_context_device context) {

	if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT_DEVICE) {
		return false;
	}

	qmckl_context_struct_device *const ctx = (qmckl_context_struct_device *)context;
	assert(ctx != NULL);

	return ctx->nucleus.provided;
}
