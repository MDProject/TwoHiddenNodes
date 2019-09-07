#include "Integrator.h"

double INT4Gaussian(double(*func)(double, double, double, double, double*), double* param, int N) {
	double t, x, u, u_prime;
	double sum = 0.;
	for (int i = 0; i < N; i++) {
		t = GaussianDistribution();
		x = GaussianDistribution();
		u = GaussianDistribution();
		u_prime = GaussianDistribution();
		sum += func(t, x, u, u_prime, param);
	}
	sum /= N;
	return sum;
}

double INT3Gaussian(double(*func)(double, double, double), int N) {
	double z1, z2, z3;
	double sum = 0.;
	for (int i = 0; i < N; i++) {
		z1 = GaussianDistribution();
		z2 = GaussianDistribution();
		z3 = GaussianDistribution();
		sum += func(z1, z2, z3);
	}
	sum /= N;
	return sum;
}

double INT2Gaussian(double(*func)(double, double, double*), double* param_hat, int N) {
	double x, y;
	double sum = 0.;
	for (int i = 0; i < N; i++) {
		x = GaussianDistribution();
		y = GaussianDistribution();
		sum += func(x, y, param_hat);
	}
	sum /= N;
	return sum;
}