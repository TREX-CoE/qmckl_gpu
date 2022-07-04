
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
CFLAGS="-Iinclude -DHAVE_HPC ${CFLAGS}"


## Enable GPU offloading

# OpenMP offloading
AC_ARG_ENABLE(openmp, [AS_HELP_STRING([--enable-openmp],[Enable OpenMP offloaded functions])], HAVE_OPENMP_OFFLOAD=yes, HAVE_OPENMP_OFFLOAD=no)
AS_IF([test "$HAVE_OPENMP_OFFLOAD" = "yes"], [
  AC_DEFINE([HAVE_OPENMP_OFFLOAD], [1], [If defined, activate OpenMP-offloaded routines])

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


# Device pointers (if OpenMP is enabled)
AC_ARG_WITH(device-pointers, [AS_HELP_STRING([--with-device-pointers],[Enable OpenMP functions using device pointers to optimize memory transfers])], DEVICE_POINTERS_TMP=yes, DEVICE_POINTERS_TMP=no)
HAVE_DEVICE_POINTERS="no"
AS_IF([test "DEVICE_POINTERS_TMP" = "yes" && test "HAVE_OPENMP_OFFLOAD" = "yes"], [
  HAVE_DEVICE_POINTERS="yes"
  AC_DEFINE([HAVE_DEVICE_POINTERS], [1], [If defined, use device pointers (requires to enable OpenMP)])
])

# cuBLAS versions (if OpenMP is enabled)
AC_ARG_WITH(cublas, [AS_HELP_STRING([--with-cublas],[Enable OpenMP functions using cuBLAS ccalls])], HAVE_CUBLAS=yes, HAVE_CUBLAS=no)
# TODO


# OpenACC offloading
AC_ARG_ENABLE(openacc, [AS_HELP_STRING([--enable-openacc],[Enable OpenACC offloaded functions])], HAVE_OPENACC_OFFLOAD=yes, HAVE_OPENACC_OFFLOAD=no)
AS_IF([test "$HAVE_OPENACC_OFFLOAD" = "yes"], [
  AC_DEFINE([HAVE_OPENACC_OFFLOAD], [1], [If defined, activate OpenACC-offloaded routines])

  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -fopenacc"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=acc"
          ;;

    *)
          echo "Warning : compiler not supported for OpenMP offload. Make sure to define necessary flags manually."
  esac
])



# Fallback solution if nothing is enabled :
# (probably because no GPU argument have been provided),
# Enable OpenMP by default
AS_IF([test "$HAVE_OPENMP_OFFLOAD" = "no"], [
  AS_IF([test "$HAVE_OPENACC_OFFLOAD" = "no"], [

    AC_OUTPUT
    echo \
    "Warning: No GPU functions could be enabled using the passed arguments. Enabling OpenMP offload by default."


    HAVE_OPENMP_OFFLOAD="yes"
    AC_DEFINE([HAVE_OPENMP_OFFLOAD], [1], [If defined, activate OpenMP-offloaded routines])

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

OpenMP offload.....: ${HAVE_OPENMP_OFFLOAD}
  Device pointers..: ${HAVE_DEVICE_POINTERS}
  cuBLAS versions..: ${HAVE_CUBLAS}
OpenACC offload....: ${HAVE_OPENACC_OFFLOAD}


Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build ${PACKAGE_NAME}
  TODO check    - run tests
  install       - install
--------------------------------------------------"

# TODO Write tests? But this requires linking with an already installed QMCkl,
# which is not necessary when simply building the library
