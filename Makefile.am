###
# Build sources
###


ACLOCAL_AMFLAGS = -I m4 -I include 

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include -I$(top_srcdir)/qmckl/src -I$(top_srcdir)/qmckl/include

include_HEADERS = include/qmckl_gpu.h

if HAVE_FORTRAN_INTERFACE
  include_HEADERS += include/qmckl_gpu_f.f90
endif

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

# (Legacy) Enable OpenMP offload
if HAVE_OMP_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_omp_offload.c src/qmckl_mo_omp_offload.c 
endif

# (Legacy) Enable OpenACC offload
if HAVE_ACC_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_acc_offload.c src/qmckl_mo_acc_offload.c 
endif


# Base library files (needed for both OpenMP and OpenACC)
if HAVE_DEVICE_POINTERS
libqmckl_gpu_la_SOURCES += src/qmckl_context.c src/qmckl_memory.c src/qmckl_blas.c \
src/qmckl_trexio.c src/qmckl_ao.c src/qmckl_mo.c
endif

# Enable OpenMP-specific device functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_omp.c \
  src/qmckl_blas_omp.c src/qmckl_trexio_omp.c \
  src/qmckl_ao_omp.c src/qmckl_mo_omp.c
endif

# Enable OpenACC-specific device functions
if HAVE_ACC_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_acc.c \
  src/qmckl_blas_acc.c src/qmckl_trexio_acc.c \
  src/qmckl_ao_acc.c src/qmckl_mo_acc.c
endif

# Enable Fortran interface
if HAVE_FORTRAN_INTERFACE
  libqmckl_gpu_la_SOURCES += src/qmckl_gpu_f.c
endif


###
# Build test programs
###

check_PROGRAMS = tests/test_qmckl_ao tests/test_qmckl_mo

qmckl/src/libqmckl.la:
	cd qmckl && make && cd ..

if HAVE_OMP_DEVICE
  tests_test_qmckl_ao_SOURCES = tests/test_qmckl_ao_omp.c chbrclf.h
  tests_test_qmckl_ao_LDFLAGS = -lm

  tests_test_qmckl_mo_SOURCES = tests/test_qmckl_mo_omp.c chbrclf.h
  tests_test_qmckl_mo_LDFLAGS = -lm
endif

if HAVE_ACC_DEVICE
  tests_test_qmckl_ao_SOURCES = tests/test_qmckl_ao_acc.c chbrclf.h
  tests_test_qmckl_ao_LDFLAGS = -lm

  tests_test_qmckl_mo_SOURCES = tests/test_qmckl_mo_acc.c chbrclf.h
  tests_test_qmckl_mo_LDFLAGS = -lm
endif

# TODO
#test_qmckl_mo_SOURCES = tests/test_qmckl_mo.c
#test_qmckl_mo_LDFLAGS = -lqmckl -lqmckl_gpu

# This shoud be added whether we use OpenMP or ACC
tests_test_qmckl_ao_LDADD = qmckl/src/libqmckl.la libqmckl_gpu.la
tests_test_qmckl_mo_LDADD = qmckl/src/libqmckl.la libqmckl_gpu.la


###
# Run test programs
###

TESTS = $(check_PROGRAMS)
