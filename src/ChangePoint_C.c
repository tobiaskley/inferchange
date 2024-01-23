#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(intervalF)(int *n2, double *dec, int *dep2, double *ilen, int *nint2, int *nsum2, int *bou);

extern SEXP intervalC(SEXP n, SEXP dec, SEXP dep, SEXP ilen, SEXP nint, SEXP nsum){
  SEXP bou;
  int *n2 = INTEGER(n);
  int *dep2 = INTEGER(dep);
  int *nint2 = INTEGER(nint);
  int *nsum2 = INTEGER(nsum);
  PROTECT(bou = allocMatrix(INTSXP, *nsum2, 2));
  F77_CALL(intervalF)(n2, REAL(dec), dep2, REAL(ilen), nint2, nsum2, INTEGER(bou));
  UNPROTECT(1);
  return(bou);
}

static const R_CallMethodDef CallEntries[] = {
  {"intervalC",   (DL_FUNC) &intervalC,   6},
  {NULL, NULL, 0}
};


void R_init_ChangePoints(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


