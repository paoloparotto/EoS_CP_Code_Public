#ifndef ROOT_FUNC_H
#define ROOT_FUNC_H

/* declare parameters and functions */
struct snB_params{
      double a;
};


double splined_func (double x, void *params);

double splined_func_deriv (double x, void *params);

void spline_func_fdf (double x, void *params, double *y, double *dy);

#endif