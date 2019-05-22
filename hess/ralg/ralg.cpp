#include "ralg.h"
#include <cfloat>
#include <mkl_cblas.h>
#include <malloc.h>
#include <cstdio>
#include <ctime>

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

double** dalloc(unsigned int n)
{
  // TODO casts
  double **all_x = (double**) malloc(n * sizeof *all_x );
  if(all_x == NULL)
    printf("allocation failed (%u)\n", n);
  double  *all_y = (double*)  calloc(n * n,  sizeof *all_y );
  if(all_y == NULL)
    printf("allocation_failed (%u)\n", n);
  double **result = all_x;
  unsigned int i;
  for (i = 0 ; i < n ; i++, all_y += n )
    result[i] = all_y;
  return result;
}

void dfree(double** m)
{
  free(m[0]);
  free(m);
}

double ralg(const ralg_options* opt,
          std::function<bool (const double*, double&, double*)> cb_grad_and_func,
          unsigned int DIMENSION,
          double* x0,
          double* res,
          bool is_min)
{
  double* xk;
  double** B;
  double* grad;
  double* tmp; // used for different tasks, store one for memory reduce
  double* tmp2;

  unsigned int i,j;
  unsigned int iter = 0;
  double step = opt->initstep;
  double d_var;
  double step_diff;
  double f_val;

  double f_optimal;

  unsigned int nr_matrix_reset = 0;
  printf("Running ralg_blas v2 with matrix renewal, copyright Eugene Lykhovyd, 2014-2018.\n");

  time_t t_started = time(NULL);

  if(opt->b_init <= 0)
  {
    printf("opt->b_init wrong value %e\n", opt->b_init);
    return 0.;
  }
  // null after init
  B = dalloc(DIMENSION);
  for(i = 0; i < DIMENSION; ++i)
    B[i][i] = opt->b_init*1.;

  xk = (double*) malloc(sizeof(double)*DIMENSION);
  grad = (double*) malloc(sizeof(double)*DIMENSION);
  tmp = (double*) malloc(sizeof(double)*DIMENSION);
  tmp2 = (double*) malloc(sizeof(double)*DIMENSION);

  cblas_dcopy(DIMENSION, x0, 1, xk, 1);
  printf("init done\n");
  time_t t_inited = time(NULL);
  if(!cb_grad_and_func(xk, f_val, grad))
  {
    printf("grad failed, aborting\n");
    return 0.;
  }

  f_optimal = f_val;

  do
  {
    iter++;

    cblas_dgemv(CblasRowMajor, CblasTrans, DIMENSION, DIMENSION, ((is_min)?(1.):(-1.)), B[0], DIMENSION, grad, 1, 0., tmp, 1);
    d_var = cblas_dnrm2(DIMENSION, tmp, 1);

    if(d_var < opt->b_mult_grad_min)
    {
      printf("B(grad) is 0, break\n");
      break;
    }

    cblas_dgemv(CblasRowMajor, CblasNoTrans, DIMENSION, DIMENSION, 1./d_var, B[0], DIMENSION, tmp, 1, 0., tmp2, 1);
    // now tmp2 is the vector we are moving in direction to
    // running adaprive step
    i=0;
    j=0;

    d_var = cblas_dnrm2(DIMENSION, tmp2, 1);
    // save previous gradient
    cblas_dcopy(DIMENSION, grad, 1, tmp, 1);

    step_diff = 0.;

    do
    {
      // tmp2 - min direction
      // tmp - old gradient
      // grad - new graient
      i++; j++;
      cblas_daxpy(DIMENSION, -step, tmp2, 1, xk, 1);

      step_diff = step_diff + step;

      if(!cb_grad_and_func(xk, f_val, grad))
      {
        printf("grad failed\n");
        break;
      }
      if(i == opt->nh)
      {
        step = step * opt->q2;
        i = 0;
      }
      if(((is_min)?(1.):(-1.))*cblas_ddot(DIMENSION, grad, 1, tmp2, 1) <= 0.)
          break;
      if(j > opt->stepmax)
      {
        printf("function is unbounded, done %d steps, current step %.14e\n", j, step);
        return 0.;
      }
    } while(1);

    if(!opt->is_monotone)
    {
      if((is_min && f_val < f_optimal) || (!is_min && f_val > f_optimal))
        cblas_dcopy(DIMENSION, xk, 1, res, 1);
    }

    if(is_min)
      f_optimal = min(f_optimal, f_val);
    else
      f_optimal = max(f_optimal, f_val);

    step_diff = step_diff * d_var;
    if(step_diff < opt->stepmin)
    {
      printf("step is 0, break\n");
      break;
    }

    if(j == 1)
      step = step * opt->q1; //decreasing

    cblas_daxpy(DIMENSION, -1., grad, 1, tmp, 1);
    cblas_dgemv(CblasRowMajor, CblasTrans, DIMENSION, DIMENSION, ((is_min)?(-1.):(1)), B[0], DIMENSION, tmp, 1, 0., tmp2, 1);
    d_var = cblas_dnrm2(DIMENSION, tmp2, 1);
    if (opt->output && (iter-1) % opt->output_iter == 0)
    {
        printf("iter = %d, step = %.14e, func = %.14e, norm = %.14e, diff = %.14e\n", iter, step, f_val, d_var, step_diff);
    }
    if(d_var > opt->reset)
    {
      cblas_dscal(DIMENSION, 1./d_var, tmp2, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, DIMENSION, DIMENSION, 1., B[0], DIMENSION, tmp2, 1, 0., tmp, 1);
      cblas_dger(CblasRowMajor, DIMENSION, DIMENSION, (1. / opt->alpha - 1.), tmp, 1, tmp2, 1, B[0], DIMENSION);
    }
    else
    {
      printf("Matrix reset on iter %d\n", iter);

      nr_matrix_reset ++;
      cblas_dscal(DIMENSION*DIMENSION, 0, B[0], 1);
      for(i=0;i<DIMENSION;++i)
        B[i][i] = 1.;
      step = step_diff / opt->nh;
    }

    if(iter > opt->itermax)
    {
      printf("max_iter reached\n");
      break;
    }
  } while(step > opt->stepmin);
  if(step <= opt->stepmin)
  {
      printf("stepmin reached\n");
  }

  if(opt->is_monotone)
    cblas_dcopy(DIMENSION, xk, 1, res, 1);

  time_t t_done = time(NULL);

  printf("ralg done, iterations : %d, matrix resets : %d\n", iter, nr_matrix_reset);
  printf("f_optimal = %e\n", f_optimal);
  printf("Time stats : init %.1lf, compute %.1lf, total %.1lf\n", difftime(t_inited, t_started), difftime(t_done, t_inited), difftime(t_done, t_started));

  //memory release
  free(tmp2);
  free(tmp);
  free(grad);
  free(xk);
  dfree(B);

  return f_optimal;
}

