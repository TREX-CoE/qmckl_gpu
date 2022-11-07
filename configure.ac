
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# QMCkl - Quantum Monte Carlo kernel library (GPU addon)
#
# BSD 3-Clause License
#
# Copyright (c) 2020, TREX Center of Excellence
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

AC_PREREQ([2.69])

AC_INIT([qmckl_gpu],[0.0.1],[https://github.com/TREX-CoE/qmckl/issues],[],[https://trex-coe.github.io/qmckl/index.html]	)
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AM_PROG_AR

AC_PROG_CC
AC_CONFIG_FILES([Makefile])
AC_OUTPUT


# Initialize Makefile
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
AC_ENABLE_SHARED
AC_ENABLE_STATIC

# Add mandatory CFLAGS
CFLAGS="-DHAVE_HPC -Iinclude ${CFLAGS}"


## Enable the different function types (OpenMP, OpenACC, device)

# OpenMP offloading
AC_ARG_ENABLE(omp_offload, [AS_HELP_STRING([--enable-omp-offload],[Enable OpenMP offloaded functions])], have_omp_offload=yes, have_omp_offload=no)

# OpenACC offloading
AC_ARG_ENABLE(acc_offload, [AS_HELP_STRING([--enable-acc-offload],[Enable OpenACC offloaded functions])], have_acc_offload=yes, have_acc_offload=no)

# OpenMP + device pointers
AC_ARG_ENABLE(omp_device, [AS_HELP_STRING([--enable-omp-device],[Enable OpenMP device pointers functions])], have_omp_device=yes, have_omp_device=no)

# OpenACC + device pointers
AC_ARG_ENABLE(acc_device, [AS_HELP_STRING([--enable-acc-device],[Enable OpenACC device pointers functions])], have_acc_device=yes, have_acc_device=no)


# TODO Add cuBLAS functions at some point ?

# Fallback solution if nothing is enabled :
# (probably because no GPU argument have been provided),
# Enable OpenMP by default
AS_IF([ test "x$have_omp_offload"  = "xno" -a "x$have_acc_offload" = "xno" -a  "x$have_omp_device"  = "xno" -a "x$have_acc_device"  = "xno"], [
	AC_OUTPUT 
    echo \
    "Warning: No GPU functions could be enabled using the passed arguments. Enabling OpenMP offload by default."

    have_omp_offload="yes"
])

AM_CONDITIONAL([HAVE_OMP_OFFLOAD], [test "x$have_omp_offload"  = "xyes"])
AM_CONDITIONAL([HAVE_ACC_OFFLOAD], [test "x$have_acc_offload" = "xyes"])
AM_CONDITIONAL([HAVE_OMP_DEVICE],  [test "x$have_omp_device"  = "xyes"])
AM_CONDITIONAL([HAVE_ACC_DEVICE],  [test "x$have_acc_device"  = "xyes"])

## Add compilers flags depending on what has been enabled
AS_IF([ test "x$have_omp_offload" = "xyes" -o "x$have_omp_device" = "xyes" ], [
  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -fopenmp"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=gpu"
          ;;

    *)
          echo "Warning : compiler not recognized for OpenMP. Make sure to define necessary flags manually."
  esac
])

AS_IF([test "x$have_acc_offload" = "xyes" -o "x$have_acc_device" = "xyes"], [
  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -fopenacc"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=acc"
          ;;

    *)
          echo "Warning : compiler not supported for OpenACC. Make sure to define necessary flags manually."
  esac
])


# Print configure status
AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'

CC.................: ${CC}
CFLAGS.............: ${CFLAGS}
LDFLAGS:...........: ${LDFLAGS}

OpenMP offload.....: ${have_omp_offload}
OpenACC offload....: ${have_acc_offload}
OpenMP device......: ${have_omp_device}
OpenACC device.....: ${have_acc_device}


Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build ${PACKAGE_NAME}
  TODO check    - run tests
  install       - install
--------------------------------------------------"

# TODO Write tests?
