ACLOCAL_AMFLAGS = -I m4 -I include 

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include -I$(top_srcdir)/qmckl/src -I$(top_srcdir)/qmckl/include

include_HEADERS = include/qmckl_gpu.h include/qmckl_gpu_f.F90

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

# Tests
TESTS = test_qmckl_ao_f test_qmckl_mo 
check_PROGRAMS = $(TESTS)
test_qmckl_ao_f_SOURCES = tests/test_qmckl_ao.c tests/qmckl_f.F90 tests/chbrclf.h tests/test_qmckl_ao_f.F90  
test_qmckl_ao_f_LDFLAGS = -lqmckl
test_qmckl_mo_SOURCES = tests/test_qmckl_mo.c  
test_qmckl_mo_LDFLAGS = -lqmckl

# (Legacy) Enable OpenMP offload
if HAVE_OMP_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_omp_offload.c
                             # src/qmckl_jastrow_omp_offload.c
endif

# (Legacy) Enable OpenACC offload
if HAVE_ACC_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_acc_offload.c src/qmckl_jastrow_acc_offload.c
endif


# Base library files (needed for both OpenMP and OpenACC)
libqmckl_gpu_la_SOURCES += src/qmckl_context.c src/qmckl_memory.c src/qmckl_blas.c src/qmckl_trexio.c

# Enable OpenMP-specific device functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_memory_omp.c \
  src/qmckl_blas_omp.c src/qmckl_trexio_omp.c \
  src/qmckl_ao_omp.c src/qmckl_jastrow_omp.c
endif

# TODO Enable OpenACC-specific device functions
