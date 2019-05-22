#ifndef RALG_H
#define RALG_H

#include <cfloat>
#include <functional> // C++11 std::function
#include <climits>

#define RALG_MAX false
#define RALG_MIN true

#define RALG_UNLIMITED_ITER UINT_MAX

struct ralg_options
{
    double q1; // 0.95
    double q2; // 1.2
    unsigned int nh; // 3
    double alpha; // 2
    double renew_limit; // DBL_EPSILON
    unsigned int itermax; // 10000
    unsigned int stepmax; // 500
    double stepmin; // 1.e-15
    double b_mult_grad_min; // 1.e-20
    double reset; // 1.e-18
    double initstep; // 1.
    bool output;
    unsigned int output_iter;
    double b_init;
    bool is_monotone;
};

const ralg_options defaultOptions = {
  // q1
  0.95,
    // q2
  1.2,
    // nh
  3,
    // alpha
  2.,
  // renew_limit
  2*DBL_EPSILON, // use with caution
    // itermax
  10000,
    // stepmax
  500,
    // stepmin
  1.e-12,
    // b_mult_grad_min
  1.e-12,
  // reset
  1.e-12,
  // initstep
  1., // so x0 = [1, 1, 1, ... , 1]
  // output
  true,
  // output_iter
  50,
  // b_init
  1.,
  // is_monotone
  true
};

double ralg(const ralg_options* opt,
          std::function<bool (const double*, double&, double*)> cb_grad_and_func,
          unsigned int DIMENSION,
          double* x0,
          double* res, bool min=RALG_MIN);

#endif // RALG_H

