#include "root_func.h"




/* define parameters and functions */
struct snB_params{
      double a;
};

double splined_func (double x, void *params)
{
  struct snB_params *p
    = (struct snB_params *) params;

  double a = p->a;

  double func = gsl_spline2d_eval(spline, xi, yi, xacc, yacc);

  return func - a;
}

double
splined_func_deriv (double x, void *params)
{
  struct snB_params *p
    = (struct snB_params *) params;

  double a = p->a;

  double func_deriv = gsl_spline2d_eval_deriv_x(spline, xi, yi, xacc, yacc);

  return func_deriv;
}

void splined_func_fdf (double x, void *params,
               double *y, double *dy)
{
  struct snB_params *p
    = (struct snB_params *) params;

  double a = p->a;

  double func = gsl_spline2d_eval(spline, xi, yi, xacc, yacc);
  double func_deriv = gsl_spline2d_eval_deriv_x(spline, xi, yi, xacc, yacc);

  *y = func;
  *dy = func_deriv;
}

