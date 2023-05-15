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


# Add mandatory CFLAGS/FCFLAGS
CFLAGS="-Iinclude ${CFLAGS}"
FCFLAGS="-Iinclude ${FCFLAGS}"


###
# TREXIO
###

AC_ARG_ENABLE([trexio],
            [AS_HELP_STRING([--disable-trexio],
                            [disable support for TREXIO])],
            [have_trexio=$enableval],
            [have_trexio=yes])

AS_CASE([$have_trexio],
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
  [AS_IF([test ! -d "$enableval"],
         [AC_MSG_ERROR([--with-trexio=$enableval: argument is not a directory])])
   TREXIO_LIBS="-L$with_trexio/lib -ltrexio"
   TREXIO_CFLAGS="-I$with_trexio/include"
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_CHECK_HEADERS([trexio.h], [],
                    [AC_MSG_ERROR("trexio.h not found")])
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"
  ])

###
# cuBLAS / cuSolver (optional, enables additional kernels)
# TODO We should probably add compilation flags to link cuBLAS / solver automatically
# in case the flag is enabled
###

AC_ARG_ENABLE([cublas],
            [AS_HELP_STRING([--enable-cublas],
                            [Enable compilation of kernels that use cuBLAS and cuSolver])],
            [have_cublas=yes],
            [have_cublas=no])

AC_ARG_ENABLE([floatmos],
            [AS_HELP_STRING([--enable-floatmos],
                            [Use single precision to compute MOs (with --enable-cublas only)])],
            [have_floatmos=yes],
            [have_floatmos=no])

AC_ARG_ENABLE([cusparse],
            [AS_HELP_STRING([--enable-cusparse],
                            [Enable compilation of kernels that use cuSPARSE])],
            [have_cusparse=yes],
            [have_cusparse=no])

AS_IF([ test "x$have_cublas" = "xyes" -a "x$have_cusparse" = "xyes" ], [
	AC_OUTPUT
	echo \
	"Warning: cannot activate cuBLAS and cuSPARSE at the same time, activating cuBLAS only."
	have_cusparse="no"
])

AS_IF([ test "x$have_floatmos" = "xyes" -a "x$have_cublas" = "xno" ], [
	AC_OUTPUT
	echo \
	"Warning: --enable-floatmos working only with --enable-cublas, deactivating. " 
	have_floatmos="no"
])


####
## Enable OpenMP / OpenACC
####

# OpenMP + device pointers
AC_ARG_ENABLE(omp, [AS_HELP_STRING([--enable-omp],[Enable OpenMP device pointers functions])], have_omp=yes, have_omp=no)

# OpenACC + device pointers
AC_ARG_ENABLE(acc, [AS_HELP_STRING([--enable-acc],[Enable OpenACC device pointers functions])], have_acc=yes, have_acc=no)

# Auto offloading flags
AC_ARG_ENABLE(autoflags, [AS_HELP_STRING([--disable-autoflags],[By default, configure will try to infer and add offloading compiler flags by itself. Disable this option and set your flags manually if you encounter issues.])], have_autoflags=no, have_autoflags=yes)
AS_IF([test "x$have_autoflags" = "xyes"], [
	echo "Autoflags have been kept enabled"
], [
	echo "Autoflags have been disabled"
])


####
## Enable the optional Fortran interface
####

AC_ARG_ENABLE(fortran-interface, [AS_HELP_STRING([--enable-fortran-interface],[Enable the optional Fortran interface])], have_fortran_interface=yes, have_fortran_interface=no)


####
## Perform checks on what has been enabled (at least one function type, no OpenACC and OpenMP at the same time)
####

# Fallback solution if nothing is enabled :
# (probably because no GPU argument have been provided)
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
AS_IF([test "x$have_autoflags" = "xyes" -a "x$have_omp" = "xyes"], [
	case $CC in
		*gcc*)
			CFLAGS="$CFLAGS -fopenmp -foffload=-lm"
			;;
		*nvc*)
			CFLAGS="$CFLAGS -mp=gpu,noautopar"
			;;
		*clang*)
			CFLAGS="$CFLAGS -fopenmp "
			;;


		*)
			echo "Warning : compiler not recognized for OpenMP. Make sure to add necessary offload flags manually."
	esac
])

AS_IF([test "x$have_autoflags" = "xyes" -a "x$have_acc" = "xyes"], [
	case $CC in
		*gcc*)
			CFLAGS="$CFLAGS -fopenacc -foffload=-lm"
			LDFLAGS="$LDFLAGS -lm -fopenacc"
			;;
		*nvc*)
			CFLAGS="$CFLAGS -acc=gpu,noautopar -gpu=nordc"
		;;

		*)
			echo "Warning : compiler not supported for OpenACC. Make sure to add necessary offload flags manually."
	esac
])

# If TREXIO or cuBLAS have been enabled, define Automake & preprocessors variables
AS_IF([test "x$have_trexio" = "xyes"], [
	AC_DEFINE([HAVE_TREXIO], [1], [Enable compilation of TREXIO functions])
])
AS_IF([test "x$have_cublas" = "xyes"], [
	AC_DEFINE([HAVE_CUBLAS], [1], [Enable compilation of cuBLAS/cuSolver functions])
	case $CC in
		*nvc*)
			CFLAGS="$CFLAGS -cudalib=cublas,cusolver"
		;;

		*)
			echo "Warning : compiler not supported for cuBLAs linking"
	esac
])

AS_IF([test "x$have_cusparse" = "xyes"], [
	AC_DEFINE([HAVE_CUSPARSE], [1], [Enable compilation of cuSPARSE functions])
	case $CC in
		*nvc*)
			echo "Warning : cuSPARSE needs cuda>=11.8. Setting -gpu=cuda11.8, change if another version is present inthe system." 
			CFLAGS="$CFLAGS -gpu=cuda11.8 -cudalib=cusparse"
		;;

		*)
			echo "Warning : compiler not supported for cuBLAs linking"
	esac
])

AS_IF([test "x$have_floatmos" = "xyes"], [
	AC_DEFINE([HAVE_FLOATMOS], [1], [Using single-precision cuBLAS])
])

# Specific options required with some compilers
AS_CASE([$FC],
 [*ifort*], [FCFLAGS="$FCFLAGS -nofor-main"],
 [*nvfortran*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf90*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf77*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [])


####
## Define Makefile variables from configure flags
####

AM_CONDITIONAL([HAVE_OMP_DEVICE],  [test "x$have_omp"  = "xyes"])
AM_CONDITIONAL([HAVE_ACC_DEVICE],  [test "x$have_acc"  = "xyes"])
AM_CONDITIONAL([HAVE_FORTRAN_INTERFACE],  [test "x$have_fortran_interface"  = "xyes"])
AM_CONDITIONAL([HAVE_DEVICE_POINTERS], [test "x$have_omp"  = "xyes" -o "x$have_acc"  = "xyes"])
AM_CONDITIONAL([HAVE_CUBLAS], [test "x$have_cublas"  = "xyes"])

####
# Print configure status
####

AC_OUTPUT

echo \
"-------------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'

CC.................: ${CC}
FC.................: ${FC}
CFLAGS.............: ${CFLAGS}
LDFLAGS:...........: ${LDFLAGS}

Autoflags .........: ${have_autoflags}
Fortran interface .: ${have_fortran_interface}

TREXIO ............: ${have_trexio}
cuBLAS ............: ${have_cublas}
cuSPARSE ..........: ${have_cusparse}

OpenMP ............: ${have_omp}
OpenACC ...........: ${have_acc}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build ${PACKAGE_NAME}
  check         - run tests
  install       - install
--------------------------------------------------"
