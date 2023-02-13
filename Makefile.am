ACLOCAL_AMFLAGS = -I m4 -I include 

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include -I$(top_srcdir)/qmckl/src -I$(top_srcdir)/qmckl/include

include_HEADERS = include/qmckl_gpu.h include/qmckl_gpu_f.F90

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

# (Legacy) Enable OpenMP offload
if HAVE_OMP_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_omp_offload.c
endif

# (Legacy) Enable OpenACC offload
if HAVE_ACC_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_acc_offload.c src/qmckl_jastrow_acc_offload.c
endif


# Base library files (needed for both OpenMP and OpenACC)
libqmckl_gpu_la_SOURCES += src/qmckl_context.c src/qmckl_memory.c src/qmckl_blas.c \
src/qmckl_trexio.c src/qmckl_ao.c src/qmckl_mo.c

# Enable OpenMP-specific device functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_omp.c \
  src/qmckl_blas_omp.c src/qmckl_trexio_omp.c \
  src/qmckl_ao_omp.c src/qmckl_mo_omp.c
endif

# Enable OpenACC-specific device functions
if HAVE_ACC_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_acc.c \
  src/qmckl_blas_omp.c src/qmckl_trexio_acc.c \
  src/qmckl_ao_acc.c src/qmckl_mo_acc.c
endif


# Go to test Makefile
SUBDIRS = tests
