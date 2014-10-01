/* 
   Additional functions to do the others calculation.
   
   Gustavo H. Esteves
   01/10/14
*/


#include <R.h>
#include <Rmath.h>
/*
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
*/




/* Defining additional functions */

double stats_mean(double *v, int n) {
    
    int i;
    double x = 0;
    
    for(i=0; i<n; i++) 
	x += v[i];
    
    return x/n;
    
}


double stats_variance(double *v, int n) {
    
    double m, x = 0;
    int i;
    
    m = stats_mean(v, n);
    
    for(i=0; i<n; i++) 
	x += R_pow((v[i] - m), 2);
    
    return (double) x/(n-1);
    
}


double stats_covariance(double *x, double *y, int n) {
    
    double m1, m2, covar=0;
    int i;
    
    m1 = stats_mean(x, n);
    m2 = stats_mean(y, n);
    
    for(i=0; i<n; i++)
	covar += (x[i] - m1)*(y[i] - m2);
    
    return (double) covar/(n-1);
    
}

double stats_sd(double *v, int n) {
    
    return sqrt(stats_variance(v, n));
    
}


double stats_median_from_sorted_data(double *v, int n) {
    
    int i;
    
    if(n%2 != 0) 
	return v[(n-1)/2];
    
    i = n/2;
    
    return (v[i] + v[i-1])/2; 
    
}


void stats_ran_sample(double *new, int k, double *obs, int n) {
    
    int numAleat, i;
    
    for(i=0; i<k; i++) {
	numAleat = (int)((unif_rand())*(n-1));
	new[i] = obs[numAleat];
    }
    
}
