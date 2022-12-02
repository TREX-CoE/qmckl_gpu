ACLOCAL_AMFLAGS = -I m4 -I include 

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include -I$(top_srcdir)/qmckl/src

include_HEADERS = include/qmckl_gpu.h

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""


# Enable OpenMP offload
if HAVE_OMP_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_omp_offload.c
                             # src/qmckl_jastrow_omp_offload.c
endif

# Enable OpenACC offload
if HAVE_ACC_OFFLOAD
  libqmckl_gpu_la_SOURCES += src/qmckl_ao_acc_offload.c src/qmckl_jastrow_acc_offload.c
endif

# Enable general device functions
# (equivalent to HAVE_OMP_DEVICE || HAVE_ACC_DEVICE,
# since  this is not possible in Automake)
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_context_device.c src/qmckl_memory_device.c \
  src/qmckl_trexio_device.c
endif

if HAVE_ACC_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_context_device.c src/qmckl_memory_device.c \
  src/qmckl_trexio_device.c
endif

# Enable OpenMP-specific device functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_context_omp_device.c \
  src/qmckl_memory_omp_device.c src/qmckl_blas_omp_device.c \
  src/qmckl_trexio_omp_device.c src/qmckl_ao_omp_device.c
endif

# TODO Enable OpenACC-specific device functions
