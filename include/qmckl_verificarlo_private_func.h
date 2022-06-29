#ifndef QMCKL_VERIFICARLO_HPT
#define QMCKL_VERIFICARLO_HPT

#include <stdbool.h>


#ifdef VFC_CI
#include <vfc_probes.h>
extern vfc_probes * probes;
#endif

#ifdef VFC_CI
void qmckl_init_probes() __attribute__((constructor));
#endif

bool qmckl_probe(
    char * testName,
    char * varName,
    double value
);

bool qmckl_probe_check(
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
);

bool qmckl_probe_check_relative(
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
);

#ifdef VFC_CI
void qmckl_dump_probes() __attribute__((destructor));
#endif

bool qmckl_probe_f(
    char * testName,
    char * varName,
    double * value
);

bool qmckl_probe_check_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
);


bool qmckl_probe_check_relative_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
);

/* [[file:../org/qmckl_verificarlo.org::*End of files][End of files:2]] */
#endif
/* End of files:2 ends here */
