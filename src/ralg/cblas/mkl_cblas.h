// This file define routines used in ralg implementation.
// Impementations can be used when efficient MKL is not available.
#ifndef _MKL_CBLAS_
#define _MKL_CBLAS_

#ifdef __cplusplus
extern "C" {
#endif

enum CBLAS_ORDER { CblasRowMajor, CblasColMajor };
enum CBLAS_TRANSPOSE { CblasNoTrans, CblasTrans };

void cblas_dcopy(const int, const double*, const int, double*, const int);
void cblas_dgemv(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const int, const int, const double, const double*, const int, const double*, const int, const double, double*, const int);
double cblas_dnrm2(const int, const double*, const int);
void cblas_daxpy(const int, const double, const double*, const int, double*, const int);
double cblas_ddot(const int, const double*, const int, const double*, const int);
void cblas_dscal(const int, const double, double*, const int);
void cblas_dger(const enum CBLAS_ORDER, const int, const int, const double, const double*, const int, const double*, const int, double*, const int);

#ifdef __cplusplus
}
#endif

#endif
