#!/bin/sh

#Parse args

usage() {
	echo "Error: this script accepts exactly 2 parameters" >&2
	echo "Usage:" >&2
	echo "[path]/[to]/copy_qmckl_headers.sh [path to qmckl_gpu] [path to qmckl]" >&2
	exit 1
}

if [ "$#" -ne 2 ]; then
	usage
fi

qmckl_gpu_path=$1
qmckl_path=$2


# Check for src and include folders existence

if [ ! -d "$qmckl_gpu_/src" ]; then
	echo "Error: src folder could not be found in $qmckl_path" >&2
	echo "Are you sure you passed the correct path to the qmckl_gpu source tree ?" >&2
fi

if [ ! -d "$qmckl_gpu_path/include" ]; then
	echo "Error: include folder could not be found in $qmckl_gpu_path" >&2
	echo "Are you sure you passed the correct path to the qmckl_gpu source tree ?" >&2
fi


# Copy private headers
headers=(
	"qmckl_ao_private_func.h"
	"qmckl_ao_private_type.h"
	"qmckl_blas_private_func.h"
	"qmckl_blas_private_type.h"
	"qmckl_context_private_type.h"
	"qmckl_determinant_private_func.h"
	"qmckl_determinant_private_type.h"
	"qmckl_electron_private_func.h"
	"qmckl_electron_private_type.h"
	"qmckl_error_private_type.h"
	"qmckl_jastrow_private_func.h"
	"qmckl_jastrow_private_type.h"
	"qmckl_local_energy_private_func.h"
	"qmckl_local_energy_private_type.h"
	"qmckl_memory_private_func.h"
	"qmckl_memory_private_type.h"
	"qmckl_mo_private_func.h"
	"qmckl_mo_private_type.h"
	"qmckl_nucleus_private_func.h"
	"qmckl_nucleus_private_type.h"
	"qmckl_point_private_func.h"
	"qmckl_point_private_type.h"
	"qmckl_trexio_private_func.h"
	"qmckl_trexio_private_type.h"
)

for header in "${headers[@]}"; do
    cp $qmckl_path/src/$header $qmckl_gpu_path/include/$header
done
