#include "../include/qmckl_trexio.h"

// This file provides wrappers to standard QMCkl functions accessible with the
// _device suffix. Only includes functions independent of OpenMP/OpenACC syntax.

//**********
// ELECTRON/POINT GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_set_electron_num_device(qmckl_context_device context,
											  const int64_t up_num,
											  const int64_t down_num) {
	return qmckl_set_electron_num((qmckl_context)context, up_num, down_num);
}

//**********
// NUCLEUS GETTERS/SETTERS
//**********

qmckl_exit_code qmckl_get_nucleus_num_device(const qmckl_context_device context,
											 int64_t *const num) {
	return qmckl_get_nucleus_num((qmckl_context)context, num);
}

qmckl_exit_code qmckl_set_nucleus_num_device(qmckl_context_device context,
											 const int64_t num) {
	return qmckl_set_nucleus_num((qmckl_context)context, num);
}

//**********
// AO GETTERS/SETTERS
//**********

qmckl_exit_code
qmckl_get_ao_basis_ao_num_device(const qmckl_context_device context,
								 int64_t *const ao_num) {
	return qmckl_get_ao_basis_ao_num((qmckl_context)context, ao_num);
}
