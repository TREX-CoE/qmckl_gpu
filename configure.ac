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

AC_INIT([qmckl_gpu],[0.0.1],[https://github.com/TREX-CoE/qmckl_gpu/issues],[],[https://trex-coe.github.io/qmckl/index.html]	)
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AM_PROG_AR

AC_PROG_CC
AC_PROG_FC
AC_CONFIG_FILES([Makefile])
AC_OUTPUT


# Initialize Makefile
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
AC_ENABLE_SHARED
AC_ENABLE_STATIC


### Find QMCkl

#AC_CHECK_LIB([qmckl], [qmckl_context_create], [], [PKG_CHECK_MODULES([QMCKL], [qmckl]) ])
#AC_CHECK_HEADERS([qmckl.h], [], [AC_MSG_ERROR("qmckl.h not found")])

# TREXIO
AC_ARG_WITH([trexio],
            [AS_HELP_STRING([--without-trexio],
                            [disable support for TREXIO])],
            [with_trexio=$withval],
            [with_trexio=yes])

AS_CASE([$with_trexio],
  [no], [],
  [yes],
  [AC_CHECK_LIB([trexio], [trexio_open], [],
                [PKG_CHECK_MODULES([TREXIO], [trexio]) ])
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_CHECK_HEADERS([trexio.h], [],
                    [AC_MSG_ERROR("trexio.h not found")])
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"],
  [AS_IF([test ! -d "$withval"],
         [AC_MSG_ERROR([--with-trexio=$withval: argument is not a directory])])
   TREXIO_LIBS="-L$with_trexio/lib -ltrexio"
   TREXIO_CFLAGS="-I$with_trexio/include"
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_CHECK_HEADERS([trexio.h], [],
                    [AC_MSG_ERROR("trexio.h not found")])
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"
  ])


# Add mandatory CFLAGS
CFLAGS="-DHAVE_HPC -Iinclude ${CFLAGS}"

####
## Enable the different function types (OpenMP, OpenACC, device)
####

# OpenMP offloading (legacy, kept for now)
AC_ARG_ENABLE(omp_offload, [AS_HELP_STRING([--enable-omp-offload],[Enable OpenMP offloaded functions])], have_omp_offload=yes, have_omp_offload=no)

# OpenACC offloading (legacy, kept for now)
AC_ARG_ENABLE(acc_offload, [AS_HELP_STRING([--enable-acc-offload],[Enable OpenACC offloaded functions])], have_acc_offload=yes, have_acc_offload=no)

# OpenMP + device pointers
AC_ARG_ENABLE(omp, [AS_HELP_STRING([--enable-omp],[Enable OpenMP device pointers functions])], have_omp=yes, have_omp=no)

# OpenACC + device pointers
AC_ARG_ENABLE(acc, [AS_HELP_STRING([--enable-acc],[Enable OpenACC device pointers functions])], have_acc=yes, have_acc=no)

# TODO Add cuBLAS functions at some point ?


####
## Perform checks on what has been enabled (at least one function type, no OpenACC and OpenMP at the same time)
####

# Fallback solution if nothing is enabled :
# (probably because no GPU argument have been provided),
# Enable OpenMP offload by default
AS_IF([ test "x$have_omp" = "xno" -a "x$have_acc" = "xno"], [
	AC_OUTPUT
	echo \
	"Warning: Neither OpenMP or OpenACC could be enabled using the passed arguments. Attempting to enable OpenMP as fallback."
	have_omp="yes"
])

# Make sure that OpenMP and OpenACC are not enabled at the same time
AS_IF([test "x$have_omp" = "xyes" -a  "x$have_acc" = "xyes"], [
	AC_OUTPUT
	echo \
	"Error: Enabling both OpenMP and OpenACC at the same time is not enabled"
	exit 1
])


####
## Set compiler flags
####

## Add compilers flags depending on what has been enabled
AS_IF([test "x$have_omp_offload" = "xyes" -o "x$have_omp" = "xyes"], [
	case $CC in
		*gcc*)
			CFLAGS="$CFLAGS -fopenmp -foffload=-lm"
			;;
		*nvc*)
			CFLAGS="$CFLAGS -mp=gpu"
			;;

		*)
			echo "Warning : compiler not recognized for OpenMP. Make sure to add necessary offload flags manually."
	esac
])

AS_IF([test "x$have_acc_offload" = "xyes" -o "x$have_acc" = "xyes"], [
	case $CC in
		*gcc*)
			CFLAGS="$CFLAGS -fopenacc -foffload=-lm"
			LDFLAGS="$LDFLAGS -lm"
			;;
		*nvc*)
			CFLAGS="$CFLAGS -acc=gpu -ta=nvidia:nordc"
		;;

		*)
			echo "Warning : compiler not supported for OpenACC. Make sure to add necessary offload flags manually."
	esac
])

# Specific options required with some compilers
AS_CASE([$FC],
 [*ifort*], [FCFLAGS="$FCFLAGS -nofor-main"],
 [*nvfortran*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf90*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf77*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [])


####
## Define Makefile variables
####

AM_CONDITIONAL([HAVE_OMP_OFFLOAD], [test "x$have_omp_offload"  = "xyes"])
AM_CONDITIONAL([HAVE_ACC_OFFLOAD], [test "x$have_acc_offload" = "xyes"])
AM_CONDITIONAL([HAVE_OMP_DEVICE],  [test "x$have_omp"  = "xyes"])
AM_CONDITIONAL([HAVE_ACC_DEVICE],  [test "x$have_acc"  = "xyes"])

####
# List other Makefiles to generate
####

AC_CONFIG_FILES([
                 tests/Makefile
				])


####
# Print configure status
####

AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'

CC.................: ${CC}
CFLAGS.............: ${CFLAGS}
LDFLAGS:...........: ${LDFLAGS}

OpenMP ............: ${have_omp}
OpenACC ...........: ${have_acc}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build ${PACKAGE_NAME}
  check         - run tests
  install       - install
--------------------------------------------------"
