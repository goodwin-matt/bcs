#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bcs_FastLaplace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bcs_intersect(SEXP, SEXP);
extern SEXP bcs_setdiff(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bcs_FastLaplace", (DL_FUNC) &bcs_FastLaplace, 6},
    {"bcs_intersect",   (DL_FUNC) &bcs_intersect,   2},
    {"bcs_setdiff",     (DL_FUNC) &bcs_setdiff,     2},
    {NULL, NULL, 0}
};

void R_init_bcs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
