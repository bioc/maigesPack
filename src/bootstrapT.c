/* 
   Function to calculate p-values from t statistic by bootstrap.
   
   Gustavo H. Esteves
   22/05/2007
*/


#include <R.h>
#include "stats.h"


/* Defining a function to calculate t statistics */
double t_stat(double *x1, double *x2, double n1, double n2, int type) {
    
    double t_res, m1, m2, s1, s2, df, v, stderro;
    
    m1 = stats_mean(x1, n1);
    m2 = stats_mean(x2, n2);
    s1 = stats_variance(x1, n1);
    s2 = stats_variance(x2, n2);
    if(type) { /* equal variances */
	df = n1+n2-2;
	v = ((n1-1)*s1 + (n2-1)*s2)/df;
	stderro = sqrt(v*(1/n1+1/n2));
    }
    else /* unequal variances */
	stderro = sqrt(s1/n1 + s2/n2);

    t_res = (m1 - m2)/stderro;  
    return t_res;
    
}

/* Defining a function to do bootstrap (bootT) */
void bootT(double *x, int *k, int *obs1, int *obs2, int *ng, int *n1, int *n2,
	   int *type, double *pval) {
    
    int i, j, nt=*n1+*n2, *resP;
    double *meant, *tstat, *m_obs1, *m_obs2, **obs1B, **obs2B, *new1, *new2,
	tstatB, *tmp, *tmp1, *tmp2;
    
    resP = (int *) R_alloc((*ng), sizeof(int));
    meant = (double *) R_alloc((*ng), sizeof(double));
    tstat = (double *) R_alloc((*ng), sizeof(double));
    m_obs1 = (double *) R_alloc((*ng), sizeof(double));
    m_obs2 = (double *) R_alloc((*ng), sizeof(double));
    new1 = (double *) R_alloc((*n1), sizeof(double));
    new2 = (double *) R_alloc((*n2), sizeof(double));
    tmp = (double *) R_alloc(nt, sizeof(double));
    tmp1 = (double *) R_alloc((*n1), sizeof(double));
    tmp2 = (double *) R_alloc((*n2), sizeof(double));
    
    /* Defining the randon generators */
    GetRNGstate();
    
    /* Defining the dimensions of the arrays */
    obs1B = (double **) R_alloc((*ng), sizeof(double));
    obs2B = (double **) R_alloc((*ng), sizeof(double));
    
    for(i=0; i<*ng; i++) {
	resP[i] = 0;
	obs1B[i] = (double *) R_alloc((*n1), sizeof(double));
	obs2B[i] = (double *) R_alloc((*n2), sizeof(double));
    }
    
    for(i=0; i < *ng; i++) {
	/* calculating the original means and medians of the data. */
	for(j=0; j<nt; j++)
	    tmp[j] = x[(nt*i)+j];
	meant[i] = stats_mean(tmp, nt);
	for(j=0; j<*n1; j++)
	    tmp1[j] = tmp[obs1[j]-1];
	m_obs1[i] = stats_mean(tmp1, *n1);
	for(j=0; j<*n2; j++)
	    tmp2[j] = tmp[obs2[j]-1];
	m_obs2[i] = stats_mean(tmp2, *n2);
	/* Calculating the original T statistic and mean (median)
	   differences. */ 
	tstat[i] = fabs(t_stat(tmp1, tmp2, *n1, *n2, *type));
	/* Generating the sequence of observations to boot. */
	for(j=0; j<*n1; j++)
	    obs1B[i][j] = tmp1[j] - m_obs1[i] + meant[i];
	for(j=0; j<*n2; j++)
	    obs2B[i][j] = tmp2[j] - m_obs2[i] + meant[i];
    }
    
    for (j=0; j<*k; j++) {
	
	if(j%500 == 0)
	    Rprintf("Doing the boots from %d up to %d\n", j+1, j+500);
	
	/* Generating the bootstraped observations */
	for(i=0; i<*ng; i++) {
	    stats_ran_sample(new1, *n1, obs1B[i], *n1);
	    stats_ran_sample(new2, *n2, obs2B[i], *n2);
	    
	    tstatB = fabs(t_stat(new1, new2, *n1, *n2, *type));
	    if(tstatB >= tstat[i])
		resP[i] += 1;
	}
	
    }

    PutRNGstate();
    
    /* Calculating the proportion of unexpected results */
    for(i=0; i<*ng; i++) 
	pval[i] = (double)resP[i]/(double)(*k);
    
}
