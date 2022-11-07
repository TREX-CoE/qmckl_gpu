ACLOCAL_AMFLAGS = -I m4 -I include

include_HEADERS = include/qmckl_gpu.h

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""


# Enable OpenMP offload
if HAVE_OMP_OFFLOAD
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openmp.c src/qmckl_jastrow_openmp.c
endif

# Enable OpenACC offload
if HAVE_ACC_OFFLOAD
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openacc.c src/qmckl_jastrow_openacc.c
endif

# Enable general device functions
if HAVE_OMP_DEVICE
 libqmckl_gpu_la_SOURCES += src/qmckl_context_device.c src/qmckl_memory_device.c \
 src/qmckl_trexio_device.c
endif

# Enable OpenMP specific functions
if HAVE_OMP_DEVICE
  libqmckl_gpu_la_SOURCES += src/qmckl_context_omp_device.c \
  src/qmckl_memory_omp_device.c src/qmckl_blas_omp_device.c \
  src/qmckl_trexio_omp_device.c src/qmckl_ao_omp_device.c
endif
