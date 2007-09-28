/* 
   Function to calculate robust correlation coeficients (minimizing the
   influences of one outlier). This strategy uses an idea similar to
   leave-one-out of cross-validation, where we make a decision based on the
   "rule of three numbers" thought by Diana an me.
   
   Gustavo H. Esteves
   22/05/07 
*/


#include <R.h>
#include "stats.h"


void corR(double *x, double *y, int n, double *res) {
    
    int i, j, k, *idx;
    double *tmp, *tmp1, *tmp2, sd1, sd2, median;
    
    idx = (int *) R_alloc(n, sizeof(int));
    tmp = (double *) R_alloc(n, sizeof(double));
    tmp1 = (double *) R_alloc((n-1), sizeof(double));
    tmp2 = (double *) R_alloc((n-1), sizeof(double));
    
    for(j=0; j<n; j++) {
	k = 0;
	for(i=0; i<n; i++)
	    if(i != j) {
		tmp1[k] = x[i];
		tmp2[k] = y[i];
		k += 1;
	    }
	sd1 = stats_sd(tmp1, (n-1));
	sd2 = stats_sd(tmp2, (n-1));
	tmp[j] = stats_covariance(tmp1, tmp2, (n-1))/(sd1*sd2);

	idx[j] = j;

    }
    
    /* Making the decision based on the rule of three numbers */

    rsort_with_index(tmp, idx, n);
    
    median = stats_median_from_sorted_data(tmp, n);
    
    if((tmp[n-1]-median) < (median-tmp[0])) {
	res[0] = tmp[0];
	res[1] = (double)idx[0];
    }
    else {
	res[0] = tmp[n-1];
	res[1] = (double)idx[n-1];
    }
    
}


void robustCorr(double *x, int *ng, int *ns, double *c, int *idx) {
  
    int i, j, k;
    double tmp[2], *tmp1, *tmp2;
    
    tmp1 = (double *) R_alloc((*ns), sizeof(double));
    tmp2 = (double *) R_alloc((*ns), sizeof(double));
    
    for (i=0; i<(*ng-1); i++)
	for(j=(i+1); j<*ng; j++) {
	    for(k=0; k<*ns; k++) {
		tmp1[k] = x[((*ns)*i)+k];
		tmp2[k] = x[((*ns)*j)+k];
	    }
	    corR(tmp1, tmp2, *ns, tmp);
	    c[((*ng)*i)+j] = tmp[0];
	    c[((*ng)*j)+i] = tmp[0];
	    idx[((*ng)*i)+j] = (int)tmp[1]+1;
	    idx[((*ng)*j)+i] = (int)tmp[1]+1;
	}
    
}
