double stats_mean(double *v, int n);
double stats_variance(double *v, int n);
double stats_covariance(double *x, double *y, int n);
double stats_sd(double *v, int n);
double stats_median_from_sorted_data(double *v, int n);
void stats_ran_sample(double *new, int k, double *obs, int n);
