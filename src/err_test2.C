####include <codecogs/stats/distributions/dists.h>
#include <algebra_polynomial_polynomial.h>
#include <maths_arithmetic_exp.h>
#include <maths_special_error_fn.h>
#include <stdio.h>
using namespace Stats::Distributions;
int main(  )
{
  double x = 0.5;
  printf("\n errorFn(%f) = %f", x, errorFn(x));
  return getchar();
}