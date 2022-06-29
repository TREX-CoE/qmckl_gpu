
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# QMCkl - Quantum Monte Carlo kernel library
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
CFLAGS="${CFLAGS} -Iinclude -DHAVE_HPC"


## Enable GPU offloading

# GPU offloading
AC_ARG_ENABLE(gpu, [AS_HELP_STRING([--enable-gpu],[openmp|openacc : Determine which GPU-offloaded functions to enable])], enable_gpu=$enableval, enable_gpu=no)
AS_IF([test "$enable_gpu" = "yes"], [enable_gpu="openmp"]) # yes equals to OpenMP

# OpenMP offloading
HAVE_OPENMP_OFFLOAD="no"
AS_IF([test "$enable_gpu" = "openmp"], [
  AC_DEFINE([HAVE_OPENMP_OFFLOAD], [1], [If defined, activate OpenMP-offloaded routines])
  HAVE_OPENMP_OFFLOAD="yes"

  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -fopenmp"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=gpu"
          ;;

    *)
          echo "Warning : compiler not supported for OpenMP offload. Make sure to define necessary flags manually."
  esac
])

# OpenACC offloading
HAVE_OPENACC_OFFLOAD="no"
AS_IF([test "$enable_gpu" = "openacc"], [
  AC_DEFINE([HAVE_OPENACC_OFFLOAD], [1], [If defined, activate OpenACC-offloaded routines])
  HAVE_OPENACC_OFFLOAD="yes"

  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -fopenacc"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -acc=gpu"
          ;;

    *)
          echo "Warning : compiler not supported for OpenACC offload. Make sure to define necessary flags manually."
  esac
])

# Device pointers
# BUG Even if the AC_DEFINE should execute only when gpu=openmp and device-pointers are set,
# the device pointers will appear as defined in the output if only device_pointers is set
# Maybe this could be related to autoconf vars being cached ?
AC_ARG_WITH(device-pointers, [AS_HELP_STRING([--with-device-pointers],[Use device pointers on offloaded functions to optimize memory transfers])], HAVE_DEVICE_POINTERS=$withval, HAVE_DEVICE_POINTERS=no)
AS_IF([test "$HAVE_DEVICE_POINTERS" = "yes" && test "$enable_gpu" = "openmp"], [
  AC_DEFINE([HAVE_DEVICE_POINTERS], [1], [If defined, use device pointers (requires to enable OpenMP)])
])


# Print configure status
AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'

CC..............: ${CC}
CFLAGS..........: ${CFLAGS}
LDFLAGS:........: ${LDFLAGS}

OpenMP offload...: ${HAVE_OPENMP_OFFLOAD}
OpenACC offload..: ${HAVE_OPENACC_OFFLOAD}
Device pointers..: ${HAVE_DEVICE_POINTERS}


Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests [TODO]
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"

# TODO Write tests? But this requires linking with an already installed QMCkl,
# which is not necessary when simply building the library
