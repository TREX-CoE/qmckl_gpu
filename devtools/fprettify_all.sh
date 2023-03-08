#!/bin/sh

#Parse args

usage() {
	echo "Error: this script accepts exactly 1 parameter" >&2
	echo "Usage:" >&2
	echo "[path]/[to]/fprettify_all.sh [path to qmckl_gpu]" >&2
	exit 1
}

if [ "$#" -ne 1 ]; then
	usage
fi

qmckl_gpu_path=$1


# Check for src and include folder existence

if [ ! -d "$qmckl_gpu_path/src" ]; then
	echo "Error: src folder could not be found in $qmckl_gpu_path" >&2
	echo "Are you sure you passed the correct path to the qmckl_gpu source tree ?" >&2
fi

if [ ! -d "$qmckl_gpu_path/include" ]; then
	echo "Error: include folder could not be found in $qmckl_gpu_path" >&2
	echo "Are you sure you passed the correct path to the qmckl_gpu source tree ?" >&2
fi



# Apply formatting
fprettify -i 4 $qmckl_gpu_path/include/qmckl_gpu_f.f90
fprettify -i 4 $qmckl_gpu_path/tests/test_qmckl_fortran_interface.f90
