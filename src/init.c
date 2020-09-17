#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(firthcox)(void *, void *, void *);
extern void F77_NAME(plcomp)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"firthcox", (DL_FUNC) &F77_NAME(firthcox), 3},
    {"plcomp",   (DL_FUNC) &F77_NAME(plcomp),   3},
    {NULL, NULL, 0}
};

void R_init_coxphf(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
