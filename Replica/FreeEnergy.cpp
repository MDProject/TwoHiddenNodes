#include "FreeEnergy.h"

double CoefA(double z1, double z3, double* param_hat, double xi1_true, double xi2_true) {
	return sqrt(param_hat[q1] - param_hat[r] / 2.)*z1 + sqrt(param_hat[r] / 2.)*z3 + param_hat[T1] * xi1_true + param_hat[tau2] * xi2_true;
}

double CoefB(double z2, double z3, double* param_hat, double xi1_true, double xi2_true) {
	return sqrt(param_hat[q2] - param_hat[r] / 2.)*z2 + sqrt(param_hat[r] / 2.)*z3 + param_hat[tau1] * xi1_true + param_hat[T2] * xi2_true;
}

double CoefC(double* param_hat) {
	return param_hat[R] - param_hat[r] / 2.;
}

double Zeff(double* param_hat, double xi1_true, double xi2_true, double z1, double z2, double z3) {
	double expt_00 = CoefA(z1, z3, param_hat, xi1_true, xi2_true)*(-1) + CoefB(z2, z3, param_hat, xi1_true, xi2_true)*(-1) + CoefC(param_hat) + log(probability_priori(-1., -1., q));
	double expt_01 = CoefA(z1, z3, param_hat, xi1_true, xi2_true)*(-1) + CoefB(z2, z3, param_hat, xi1_true, xi2_true) + CoefC(param_hat)*(-1) + log(probability_priori(-1., 1., q));
	double expt_10 = CoefA(z1, z3, param_hat, xi1_true, xi2_true) + CoefB(z2, z3, param_hat, xi1_true, xi2_true)*(-1) + CoefC(param_hat)*(-1) + log(probability_priori(1., -1., q));
	double expt_11 = CoefA(z1, z3, param_hat, xi1_true, xi2_true) + CoefB(z2, z3, param_hat, xi1_true, xi2_true) + CoefC(param_hat) + log(probability_priori(1., 1., q));
	double zeff = exp(expt_00) + exp(expt_01) + exp(expt_10) + exp(expt_11);
	return zeff;
}

double LnZeff(double z1, double z2, double z3, double* param_hat) {
	double zeff_00 = Zeff(param_hat, -1., -1., z1, z2, z3);
	double zeff_01 = Zeff(param_hat, -1., 1., z1, z2, z3);
	double zeff_10 = Zeff(param_hat, 1., -1., z1, z2, z3);
	double zeff_11 = Zeff(param_hat, 1., 1., z1, z2, z3);
	double lnZeff = log(zeff_00)*probability_priori(-1., -1., q) + log(zeff_01)*probability_priori(-1., 1., q) + log(zeff_10)*probability_priori(1., -1., q) + log(zeff_11)*probability_priori(1., 1., q);
	return lnZeff;
}

double Fenerg1(double* param_hat, int N) {
	double f_energ = 0.;
	for (int n = 0; n < N; n++) {
		double z1 = GaussianDistribution();
		double z2 = GaussianDistribution();
		double z3 = GaussianDistribution();
		f_energ += LnZeff(z1, z2, z3, param_hat);
	}
	return f_energ / N;
}

double CoefL(double* param) {
	return param[T1] * param[T1] + (param[tau2] - param[T1] * q)*(param[tau2] - param[T1] * q) / (1. - q * q);
}

double CoefH(double* param) {
	return param[T1] * param[tau1] + (param[T2] - param[tau1] * q)*(param[tau2] - param[T1] * q) / (1. - q * q);
}

double CoefS(double* param) {
	return param[tau1] * param[tau1] + (param[T2] - param[tau1] * q)*(param[T2] - param[tau1] * q) / (1. - q * q);
}

double OmegaTilde(double z, double x, double omega, double omega_tilde, double* param) {
	double L = CoefL(param);
	double H = CoefH(param);
	double S = CoefS(param);
	double term1 = (param[T1] + param[tau1])*z;
	double term2 = (param[tau2] - param[T1] * q + param[T2] - param[tau1] * q) / sqrt(1. - q * q)*x;
	double term3 = (sqrt(param[q1] - L) + (param[r] - H) / sqrt(param[q1] - L))*omega;
	
	double tmp = (param[r] - H) / sqrt(param[q1] - L);
	double term4 = sqrt(param[q2] - tmp * tmp - S)*omega_tilde;
	double plus = term1 + term2 + term3 + term4;

	term1 = (param[T1] - param[tau1])*z;
	term2 = (param[tau2] - param[T1] * q - param[T2] + param[tau1] * q) / sqrt(1. - q * q)*x;
	term3 = (sqrt(param[q1] - L) - (param[r] - H) / sqrt(param[q1] - L))*omega;
	double minus = term1 + term2 + term3 - term4;

	double Omega = exp(beta*beta*(param[R] - param[r]))*cosh(beta*plus) + exp(0. - beta * beta*(param[R] - param[r]))*cosh(beta*minus);
	return Omega;
}

double Fenerg2(double* param, int N) {
	double f_energ = 0.;
	for (int n = 0; n < N; n++) {
		double z = GaussianDistribution();
		double x = GaussianDistribution();
		double omega = GaussianDistribution();
		double omega_tilde = GaussianDistribution();
		double func_integral = cosh(beta*z)*cosh(beta*(q*z + sqrt(1 - q * q)*x))*log(OmegaTilde(z, x, omega, omega_tilde, param));
		f_energ += func_integral;
	}
	f_energ /= f_energ / N;
	f_energ *= alpha / exp(beta*beta) / cosh(beta*beta*q);
	return f_energ;
}

double CalculateFreeEnergy(double* param, double* param_hat, int N) {
	double term0 = param_hat[q1] * param[q1] / 2. + param_hat[q2] * param[q2] / 2. - param_hat[q1] / 2. - param_hat[q2] / 2. - param_hat[R] * param[R] + param_hat[r] * param[r] / 2.;
	term0 = term0 - param_hat[T1] * param[T1] - param_hat[tau1] * param[tau1] - param_hat[T2] * param[T2] - param_hat[tau2] * param[tau2] + alpha * beta*beta*(1. - param[q1] / 2. - param[q2] / 2.) - alpha * log(2.*cosh(beta*beta*param[R]));

	double term1 = Fenerg1(param_hat, N);
	
	double term2 = Fenerg2(param, N);

	return term0 + term1 + term2;
}