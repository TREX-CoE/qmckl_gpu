###
# Build test programs
###


check_PROGRAMS = test_qmckl_ao test_qmckl_mo

if HAVE_OMP_DEVICE
  test_qmckl_ao_SOURCES = test_qmckl_ao_omp.c chbrclf.h
  test_qmckl_ao_LDFLAGS = -lm

  test_qmckl_mo_SOURCES = test_qmckl_mo_omp.c chbrclf.h
  test_qmckl_mo_LDFLAGS = -lm
endif

if HAVE_ACC_DEVICE
  test_qmckl_ao_SOURCES = test_qmckl_ao_acc.c chbrclf.h
  test_qmckl_ao_LDFLAGS = -lm

  test_qmckl_mo_SOURCES = test_qmckl_mo_acc.c chbrclf.h
  test_qmckl_mo_LDFLAGS = -lm
endif

# TODO
#test_qmckl_mo_SOURCES = tests/test_qmckl_mo.c
#test_qmckl_mo_LDFLAGS = -lqmckl -lqmckl_gpu

###
# All of this shoud be added whether we use OpenMP or ACC
###

test_qmckl_ao_LDADD = ../qmckl/src/libqmckl.la ../libqmckl_gpu.la
test_qmckl_mo_LDADD = ../qmckl/src/libqmckl.la ../libqmckl_gpu.la

###
# Run test programs
###

TESTS = $(check_PROGRAMS)
