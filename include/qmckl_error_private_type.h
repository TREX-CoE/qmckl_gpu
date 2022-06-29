#ifndef QMCKL_ERROR_HPT
#define QMCKL_ERROR_HPT

/* Data structure in context */

/*   The strings are declared internally with a maximum fixed size to avoid */
/*   dynamic memory allocation. */


#define  QMCKL_MAX_FUN_LEN   256
#define  QMCKL_MAX_MSG_LEN  1024

typedef struct qmckl_error_struct {

  qmckl_exit_code exit_code;
  char function[QMCKL_MAX_FUN_LEN];
  char message [QMCKL_MAX_MSG_LEN];

} qmckl_error_struct;

/* [[file:../org/qmckl_error.org::*End of files][End of files:1]] */
#endif
/* End of files:1 ends here */
