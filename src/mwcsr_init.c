#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP _mwcsr_rmwcs_solve(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mwcsr_rmwcs_solve", (DL_FUNC) &_mwcsr_rmwcs_solve, 2},
    {NULL, NULL, 0}
};

void R_init_mwcsr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

