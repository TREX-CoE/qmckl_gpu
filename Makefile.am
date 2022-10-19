ACLOCAL_AMFLAGS = -I m4 -I include

include_HEADERS = include/qmckl_gpu.h

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

if HAVE_OPENMP
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openmp.c src/qmckl_jastrow_openmp.c
endif

if HAVE_OPENACC
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openacc.c src/qmckl_jastrow_openacc.c
endif

if HAVE_DEVICE
 libqmckl_gpu_la_SOURCES += src/qmckl_context_device.c src/qmckl_memory_device.c \
 src/qmckl_blas_device.c src/qmckl_trexio_device.c src/qmckl_ao_device.c
endif

