###
# Build sources
###


ACLOCAL_AMFLAGS = -I m4 -I include 

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include 

include_HEADERS = include/qmckl_gpu.h

# Optionally distribure the Fortran interface header
if HAVE_FORTRAN_INTERFACE
  include_HEADERS += include/qmckl_gpu_f.F90
endif

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

libqmckl_gpu_la_SOURCES += src/qmckl_context.c src/qmckl_basic_functions.c src/qmckl_point.c \
src/qmckl_memory.c src/qmckl_blas.c src/qmckl_electron.c src/qmckl_nucleus.c \
src/qmckl_ao.c src/qmckl_mo.c src/qmckl_jastrow.c

# Enable TREXIO interfacing
if HAVE_TREXIO
  libqmckl_gpu_la_SOURCES += src/qmckl_trexio.c
endif

# Enable OpenMP-specific device functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_omp.c \
  src/qmckl_nucleus_omp.c src/qmckl_point_omp.c \
  src/qmckl_blas_omp.c src/qmckl_distance_omp.c \
  src/qmckl_ao_omp.c src/qmckl_mo_omp.c src/qmckl_jastrow_omp.c
endif

if HAVE_OMP_DEVICE
if HAVE_CUBLAS
    libqmckl_gpu_la_SOURCES += src/qmckl_woodbury_omp.c
endif
endif

# Enable OpenACC-specific device functions
if HAVE_ACC_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_acc.c \
  src/qmckl_nucleus_acc.c src/qmckl_point_acc.c \
  src/qmckl_blas_acc.c src/qmckl_distance_acc.c \
  src/qmckl_ao_acc.c src/qmckl_mo_acc.c src/qmckl_jastrow_acc.c
endif


###
# Build test programs
###

check_PROGRAMS = tests/test_qmckl_ao tests/test_qmckl_mo tests/test_qmckl_jastrow

if HAVE_OMP_DEVICE
  tests_test_qmckl_ao_SOURCES = tests/test_qmckl_ao_omp.c chbrclf.h
  tests_test_qmckl_ao_LDFLAGS = -lm

  tests_test_qmckl_mo_SOURCES = tests/test_qmckl_mo_omp.c chbrclf.h
  tests_test_qmckl_mo_LDFLAGS = -lm

  tests_test_qmckl_jastrow_SOURCES = tests/test_qmckl_jastrow_omp.c n2.h
  tests_test_qmckl_jastrow_LDFLAGS = -lm
endif

if HAVE_ACC_DEVICE
  tests_test_qmckl_ao_SOURCES = tests/test_qmckl_ao_acc.c chbrclf.h
  tests_test_qmckl_ao_LDFLAGS = -lm

  tests_test_qmckl_mo_SOURCES = tests/test_qmckl_mo_acc.c chbrclf.h
  tests_test_qmckl_mo_LDFLAGS = -lm

  tests_test_qmckl_jastrow_SOURCES = tests/test_qmckl_jastrow_acc.c n2.h
  tests_test_qmckl_jastrow_LDFLAGS = -lm
endif

# This shoud be added whether we use OpenMP or ACC
tests_test_qmckl_ao_LDADD = libqmckl_gpu.la
tests_test_qmckl_mo_LDADD = libqmckl_gpu.la
tests_test_qmckl_jastrow_LDADD = libqmckl_gpu.la

if HAVE_FORTRAN_INTERFACE
  check_PROGRAMS += tests/test_qmckl_fortran_interface
  tests_test_qmckl_fortran_interface_SOURCES = tests/n2_fortran_interface.c tests/test_qmckl_fortran_interface.f90 tests/chbrclf_fortran_interface.c
  tests_test_qmckl_fortran_interface_LDADD = libqmckl_gpu.la
  CLEANFILES = qmckl_gpu_f.mod
endif

###
# Run test programs
###

TESTS = $(check_PROGRAMS)
