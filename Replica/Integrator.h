#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "Random.h"


double INT4Gaussian(double(*func)(double, double, double, double, double*), double* param, int N = 1000000);
double INT3Gaussian(double(*func)(double, double, double), int N);
double INT2Gaussian(double(*func)(double, double, double*), double* param_hat, int N = 1000000);


#endif