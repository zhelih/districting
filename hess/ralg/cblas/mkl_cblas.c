#include "mkl_cblas.h"

#include <math.h>
#include <stdio.h>

void cblas_dcopy(const int __N, const double *__X, const int __incX, double *__Y, const int __incY)
{
  int xi = 0, yi = 0;
  for(; xi < __N && yi < __N; xi+=__incX, yi+=__incY)
    __Y[yi] = __X[xi];
}

void cblas_dgemv(const enum CBLAS_ORDER __Order, const enum CBLAS_TRANSPOSE __TransA, const int __M, const int __N, const double __alpha, const double *__A, const int __lda, const double *__X, const int __incX, const double __beta, double *__Y, const int __incY)
{
  if(__Order != CblasRowMajor)
  {
    fprintf(stderr, "cblas_dgemv : Order must be CblasRowMajor\n");
    return;
  }

  // y = alpha*A*x + beta*y
  int yi = 0;
  for(; yi < __N; yi += __incY)
    __Y[yi] = __beta * __Y[yi];

  //FIXME simulate actual performance of cblas
  int i, j;
  for(i = 0; i < __M; ++i)
    for(j = 0; j < __N; ++j)
    {
      double v = (__TransA == CblasNoTrans) ? __A[i*__M + j] : __A[j*__N + i];
      __Y[i] += __alpha * v * __X[j];
    }
}

double cblas_dnrm2(const int __N, const double *__X, const int __incX)
{
  double ret = 0.;
  int xi = 0;
  for(; xi < __N; xi += __incX)
    ret += __X[xi] * __X[xi];
  return sqrt(ret);
}

void cblas_daxpy(const int __N, const double __alpha, const double *__X, const int __incX, double *__Y, const int __incY)
{
  // compute y <- alpha * x + y
  int xi;
  //FIXME what does __incY do?
  for(xi = 0; xi < __N; xi+= __incX)
    __Y[xi] += __alpha * __X[xi];
}

double cblas_ddot(const int __N, const double *__X, const int __incX, const double* __Y, const int __incY)
{
  double ret = 0.;
  //FIXME again, why need incY?
  int xi;
  for(xi = 0; xi < __N; xi += __incX)
    ret += __X[xi] * __Y[xi];
  return ret;
}

void cblas_dscal(const int __N, const double __alpha, double *__X, const int __incX)
{
  int xi;
  for(xi = 0; xi < __N; xi += __incX)
    __X[xi] *= __alpha;
}

void cblas_dger(const enum CBLAS_ORDER __Order, const int __M, const int __N, const double __alpha, const double *__X, const int __incX, const double *__Y, const int __incY, double *__A, const int __lda)
{
  //FIXME simulate exactly as cblas does this
  if(__Order != CblasRowMajor)
  {
    fprintf(stderr, "cblas_dger : Order must be CblasRowMajor\n");
    return;
  }

  int i, j;
  for(i = 0; i < __M; ++i)
    for(j = 0; j < __N; ++j)
      __A[__M*i + j] += __alpha * __X[i] * __Y[j];
}
