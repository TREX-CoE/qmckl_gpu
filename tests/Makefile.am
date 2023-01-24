###
# Build test programs
###


check_PROGRAMS = test_qmckl_ao #test_qmckl_mo
test_qmckl_ao_SOURCES = test_qmckl_ao.c chbrclf.h
test_qmckl_ao_LDFLAGS = -lm -lomp
test_qmckl_ao_LDADD = ../qmckl/src/libqmckl.la ../libqmckl_gpu.la

#test_qmckl_mo_SOURCES = tests/test_qmckl_mo.c
#test_qmckl_mo_LDFLAGS = -lqmckl -lqmckl_gpu
#test_qmckl_mo_LDADD = ../src/libqmckl_gpu.la


###
# Run test programs
###

#TESTS = $(check_PROGRAMS)
