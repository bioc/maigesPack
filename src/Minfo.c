/* 
   Function to calculate Mutual information from a matrix of data.
   
   Gustavo H. Esteves
   22/05/07
*/

#include <R.h>
#include <Rmath.h>
#include "stats.h"


double dist(double *p1, double *p2) {
    
    double tmp1 = fabs(p1[0]-p2[0]);
    double tmp2 = fabs(p1[1]-p2[1]);
    
    if(tmp1 > tmp2)
	return(tmp1);
    else
	return(tmp2);
    
}

void count(double *x, double *y, int n, int point, double cut, int *result) {
    
    int i;
    
    result[0] = 0;
    result[1] = 0;
    
    for(i=0; i<n; i++) {
	if((i != point) & (fabs(x[i]-x[point]) < cut))
	    result[0] = result[0] + 1;
	if((i != point) & (fabs(y[i]-y[point]) < cut))
	    result[1] = result[1] + 1;
    }
    
}

double meanCountK(double *x, double *y, int n, int k) {
    
    double *tmp, pi[2], pj[2], *distance, result;
    int i, j, l, contagem[2];
    
    distance = (double *) R_alloc((n-1), sizeof(double));
    tmp = (double *) R_alloc(n, sizeof(double));
    
    for(i=0; i<n; i++) {
	l = 0;
	pi[0] = x[i];
	pi[1] = y[i];
	for(j=0; j<n; j++)
	    if(i != j) {
		pj[0] = x[j];
		pj[1] = y[j];
		distance[l] = dist(pi, pj);
		l = l+1;
	    }
	
	R_rsort(distance, (n-1));
	count(x, y, n, i, distance[k-1], contagem);
	
	tmp[i] = (digamma((double)(contagem[0]+1)) + 
		  digamma((double)(contagem[1]+1)))/2;
	
    }
    
    result = stats_mean(tmp, n);
    
    return(result);
    
}


void Minfo(double *x, int *ng, int *ns, int *k, double *MI) {
    
    
    double tmp, *tmp1, *tmp2;
    int i, j, l;
    
    tmp1 = (double *) R_alloc((*ns), sizeof(double));
    tmp2 = (double *) R_alloc((*ns), sizeof(double));
    
    for (i=0; i<*ng; i++)
    for(j=i; j<*ng; j++) {
	for(l=0; l<*ns; l++) {
	    tmp1[l] = x[((*ns)*i)+l];
	    tmp2[l] = x[((*ns)*j)+l];
	}
	tmp = digamma((double)*k) - meanCountK(tmp1, tmp2, *ns, *k) + 
	    digamma((double)*ns);
	MI[((*ng)*i)+j] = tmp;
	MI[((*ng)*j)+i] = tmp;
    }
    
}
