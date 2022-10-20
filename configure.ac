
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
AC_ARG_ENABLE(openmp, [AS_HELP_STRING([--enable-openmp],[Enable OpenMP offloaded functions])], have_openmp=yes, have_openmp=no)

# OpenACC offloading
AC_ARG_ENABLE(openacc, [AS_HELP_STRING([--enable-openacc],[Enable OpenACC offloaded functions])], have_openacc=yes, have_openacc=no)

# OpenMP + device pointers
AC_ARG_ENABLE(device, [AS_HELP_STRING([--enable-device],[Enable device pointers OpenMP functions])], have_device=yes, have_device=no)

# TODO Add cuBLAS functions at some point ?

# Fallback solution if nothing is enabled :
# (probably because no GPU argument have been provided),
# Enable OpenMP by default
AS_IF([ test "x$have_openmp"  = "xno" -a "x$have_openacc" = "xno" -a  "x$have_device"  = "xno" ], [ 
	AC_OUTPUT 
    echo \
    "Warning: No GPU functions could be enabled using the passed arguments. Enabling OpenMP offload by default."

    have_openmp="yes"
])

AM_CONDITIONAL([HAVE_OPENMP],  [test "x$have_openmp"  = "xyes"])
AM_CONDITIONAL([HAVE_OPENACC], [test "x$have_openacc" = "xyes"])
AM_CONDITIONAL([HAVE_DEVICE],  [test "x$have_device"  = "xyes"])

## Add compilers flags depending on what has been enabled
AS_IF([ test "x$have_openmp" = "xyes" -o "x$have_device" = "xyes" ], [
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


AS_IF([test "x$have_openacc" = "xyes"], [
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

AS_IF([ test "x$have_device" = "xyes" ], [
  case $CC in
    *gcc*)
          CFLAGS="$CFLAGS -lomp"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -lomp"
          ;;

    *)
          echo "Warning : compiler not recognized for OpenMP. Make sure to define necessary flags manually."
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

OpenMP offload.....: ${have_openmp}
OpenACC offload....: ${have_openacc}
Device offload.....: ${have_device}


Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build ${PACKAGE_NAME}
  TODO check    - run tests
  install       - install
--------------------------------------------------"

# TODO Write tests? But this requires linking with an already installed QMCkl,
# which is not necessary when simply building the library
