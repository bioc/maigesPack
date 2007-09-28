/* 
   Additional functions to do the others calculation.
   
   Gustavo H. Esteves
   20/08/07
*/


#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>



void bootT(double *x, int *k, int *obs1, int *obs2, int *ng, int *n1, int *n2,
	   int *type, double *pval);
void Minfo(double *x, int *ng, int *ns, int *k, double *MI);
void robustCorr(double *x, int *ng, int *ns, double *c, int *idx);



/* Registering functions to be used into R */
static const R_CMethodDef cMethods[] = {
    {"bootT", (DL_FUNC) &bootT, 9},
    {"Minfo", (DL_FUNC) &Minfo, 5},
    {"robustCorr", (DL_FUNC) &robustCorr, 5},
    {NULL, NULL, 0}
};


void R_init_maigesPack(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
