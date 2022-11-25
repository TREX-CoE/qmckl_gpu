// This file will contains general functions accepting a qmckl_contect_device as
// argument, hence why they need a device alternative

#include "../include/qmckl_context_device.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code qmckl_context_touch_device(const qmckl_context_device context) {
	return qmckl_context_touch((qmckl_context)context);
}

//**********
// CONTEXT CREATE
//**********
// NOTE The context destroy function frees some memory, so its implementation is
// OpenMP/OpenACC dependent

qmckl_context_device qmckl_context_create_device(int device_id) {
	qmckl_context_struct* const ctx = (qmckl_context_struct*) qmckl_context_create();

	qmckl_context_device_struct* const ds = malloc(sizeof(qmckl_context_device_struct));
	assert (ds != NULL);
	memset(ds, 0, sizeof(qmckl_context_device_struct));
       	ds->device_id = device_id;
	ctx->qmckl_extra = (void*) ds;
	return ctx;
}
