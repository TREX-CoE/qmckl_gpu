// This file will contains general functions accepting a qmckl_contect_device as
// argument, hence why they need a device alternative

#include "../include/qmckl_context_device.h"

//**********
// MISC FUNCTIONS
//**********

qmckl_exit_code
qmckl_context_touch_device(const qmckl_context_device context)
{
  return qmckl_context_touch((qmckl_context) context);
}

//**********
// CONTEXT CREATE
//**********
// NOTE The context destroy function frees some memory, so its implementation is
// OpenMP/OpenACC dependent

qmckl_context_device qmckl_context_create_device() {
  return (qmckl_context_device) qmckl_context_create();
}

