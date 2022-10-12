ACLOCAL_AMFLAGS = -I m4 -I include

include_HEADERS = include/qmckl_gpu.h include/qmckl_ao_gpu.h include/qmckl_blas_device.h include/qmckl_context_device.h include/qmckl_device_types.h include/qmckl_memory_device.h include/qmckl_trexio_device.h

lib_LTLIBRARIES = libqmckl_gpu.la

libqmckl_gpu_la_SOURCES = ""

if HAVE_OPENMP
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openmp.c
endif

if HAVE_OPENACC
 libqmckl_gpu_la_SOURCES += src/qmckl_ao_openacc.c
endif

if HAVE_DEVICE
 libqmckl_gpu_la_SOURCES += src/qmckl_context_device.c src/qmckl_memory_device.c \
 src/qmckl_blas_device.c src/qmckl_trexio_device.c src/qmckl_ao_device.c
endif

