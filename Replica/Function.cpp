#include "Function.h"
#include <math.h>
#include "Define.h"
// [param] array stands for [T1,T2,tau1,tau2,q1,q2,R,r] with no hat
// q is pre-defined in "Define.h" file as Global constant

/*
param[0] = T1
param[1] = T2
param[2] = tau1
param[3] = tau2
param[4] = q1
param[5] = q2
param[6] = R
param[7] = r
*/

double A(double* param) {
	return param[T1] * param[tau1] + (param[tau2] - param[T1] * q)*(param[T2] - param[tau1] * q) / (1 - q * q);
}

double B(double* param) {
	return param[q1] - param[T1] * param[T1] - (param[tau2] - param[T1] * q)*(param[tau2] - param[T1] * q) / (1 - q * q);
}

double C(double* param) {
	double Aval = A(param);
	return param[q2] - param[tau1] * param[tau1] - (param[T2] - param[tau1] * q)*(param[T2] - param[tau1] * q) / (1 - q * q) - (param[r] - Aval)*(param[r] - Aval) / B(param);
}

double X0(double t, double x, double u, double u_prime, double* param) {
	return t;
}
double Y0(double t, double x, double u, double u_prime, double* param) {
	return q * t + sqrt(1 - q * q)*x;
}
double X(double t, double x, double u, double u_prime, double* param) {
	double sqrt_B = sqrt(B(param));
	return param[T1] * t + (param[tau2] - param[T1] * q) / sqrt(1 - q * q)*x + sqrt_B * u;
}
double Y(double t, double x, double u, double u_prime, double* param) {
	double sqrt_C = sqrt(C(param));
	double sqrt_B = sqrt(B(param));
	return param[tau1] * t + (param[T2] - param[tau1] * q) / sqrt(1 - q * q)*x + (param[r] - A(param)) / sqrt_B * u + sqrt_C * u_prime;
}

double GsPlus(double t, double x, double u, double u_prime, double* param) {
	// (Gs1+Gs2)/(Gs3+Gs4) in order
	double Xval = X(t, x, u, u_prime, param);
	double Yval = Y(t, x, u, u_prime, param);
	double Gs1 = exp(beta*beta*(param[R] - param[r]))*sinh(beta*(Xval + Yval));
	double Gs2 = exp(-beta*beta*(param[R] - param[r]))*sinh(beta*(Xval - Yval));
	double Gs3 = exp(beta*beta*(param[R] - param[r]))*cosh(beta*(Xval + Yval));
	double Gs4 = exp(-beta*beta*(param[R] - param[r]))*cosh(beta*(Xval - Yval));
	return (Gs1 + Gs2) / (Gs3 + Gs4);
}

double GsMinus(double t, double x, double u, double u_prime, double* param) {
	double Xval = X(t, x, u, u_prime, param);
	double Yval = Y(t, x, u, u_prime, param);
	double Gs1 = exp(beta*beta*(param[R] - param[r]))*sinh(beta*(Xval + Yval));
	double Gs2 = exp(-beta * beta*(param[R] - param[r]))*sinh(beta*(Xval - Yval));
	double Gs3 = exp(beta*beta*(param[R] - param[r]))*cosh(beta*(Xval + Yval));
	double Gs4 = exp(-beta * beta*(param[R] - param[r]))*cosh(beta*(Xval - Yval));
	return (Gs1 - Gs2) / (Gs3 + Gs4);
}

double GMinus(double t, double x, double u, double u_prime, double* param) {
	double Xval = X(t, x, u, u_prime, param);
	double Yval = Y(t, x, u, u_prime, param);
	double Gs1 = exp(beta*beta*(param[R] - param[r]))*cosh(beta*(Xval + Yval));
	double Gs2 = exp(-beta * beta*(param[R] - param[r]))*cosh(beta*(Xval - Yval));
	double Gs3 = exp(beta*beta*(param[R] - param[r]))*cosh(beta*(Xval + Yval));
	double Gs4 = exp(-beta * beta*(param[R] - param[r]))*cosh(beta*(Xval - Yval));
	return (Gs1 - Gs2) / (Gs3 + Gs4);
}

double K() {
	return alpha * beta*beta*exp(-beta*beta);
}

// Accessed by Integrator
double Function_T1_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsPlusVal = GsPlus(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * sinh(beta*X0Val)*cosh(beta*Y0Val)*GsPlusVal;
}

double Function_T2_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsMinusVal = GsMinus(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * cosh(beta*X0Val)*sinh(beta*Y0Val)*GsMinusVal;
}

double Function_tau1_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsMinusVal = GsMinus(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * sinh(beta*X0Val)*cosh(beta*Y0Val)*GsMinusVal;
}

double Function_tau2_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsPlusVal = GsPlus(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * cosh(beta*X0Val)*sinh(beta*Y0Val)*GsPlusVal;
}

double Function_q1_hat(double t, double x, double u, double u_prime, double* param) {
	double GsPlusVal = GsPlus(t, x, u, u_prime, param);
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * GsPlusVal * GsPlusVal*cosh(beta*X0Val)*cosh(beta*Y0Val);
}

double Function_q2_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsMinusVal = GsMinus(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * GsMinusVal * GsMinusVal*cosh(beta*X0Val)*cosh(beta*Y0Val);
}

double Function_R_hat(double t, double x, double u, double u_prime, double* param) {
	double GMinusVal = GMinus(t, x, u, u_prime, param);
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double Kval = K();
	return Kval * (cosh(beta*X0Val)*cosh(beta*Y0Val)*GMinusVal);
}

double Function_r_hat(double t, double x, double u, double u_prime, double* param) {
	double X0Val = X0(t, x, u, u_prime, param);
	double Y0Val = Y0(t, x, u, u_prime, param);
	double GsPlusVal = GsPlus(t, x, u, u_prime, param);
	double GsMinusVal = GsMinus(t, x, u, u_prime, param);
	double Kval = K();
	return 2 * Kval*cosh(beta*X0Val)*cosh(beta*Y0Val)*GsPlusVal*GsMinusVal;
}

// [param_hat] array stands for [T1,T2,tau1,tau2,q1,q2,R,r] with hat symbol
// q is pre-defined in "Define.h" file as Global constant
double phi(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	return param_hat[r] / sqrt(param_hat[q1] * param_hat[q2]) / 2.;
}

double a1(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	return sqrt(param_hat[q1])*z1 + param_hat[T1] * xi0_1 + param_hat[tau2] * xi0_2;
}

double a2(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	double phiVal = phi(z1, z2, xi0_1, xi0_2, param_hat);
	if (phiVal > 1) {
		phiVal = 1;
	}
	return sqrt(param_hat[q2])*(phiVal*z1 + sqrt(1 - phiVal * phiVal)*z2) + param_hat[T2] * xi0_2 + param_hat[tau1] * xi0_1;
}

double a3(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	return param_hat[R] - param_hat[r] / 2.;
}

double xiGamma_1(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	double tanh_a1 = tanh(a1(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a2 = tanh(a2(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a3 = tanh(a3(z1, z2, xi0_1, xi0_2, param_hat));
	return (tanh_a1 + tanh_a2 * tanh_a3 +q*tanh_a2+q*tanh_a1*tanh_a3  ) / (1 + tanh_a1 * tanh_a2*tanh_a3 +q*tanh_a3+q*tanh_a1*tanh_a2      );
}

double xiGamma_2(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	double tanh_a1 = tanh(a1(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a2 = tanh(a2(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a3 = tanh(a3(z1, z2, xi0_1, xi0_2, param_hat));
	return (tanh_a2 + tanh_a3 * tanh_a1+q*tanh_a1+q*tanh_a2*tanh_a3   )/(1 + tanh_a1 * tanh_a2*tanh_a3  +q*tanh_a3+q*tanh_a1*tanh_a2   );
}

double xiGamma_12(double z1, double z2, double xi0_1, double xi0_2, double* param_hat) {
	double tanh_a1 = tanh(a1(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a2 = tanh(a2(z1, z2, xi0_1, xi0_2, param_hat));
	double tanh_a3 = tanh(a3(z1, z2, xi0_1, xi0_2, param_hat));
	return (tanh_a3 + tanh_a1 * tanh_a2 +q*tanh_a1*tanh_a2*tanh_a3+q   ) / (1 + tanh_a1 * tanh_a2*tanh_a3 +q*tanh_a3+q*tanh_a1*tanh_a2  );
}

// Accessed by Integrator
double Function_T1(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_1(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_1(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_1(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_1(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11 + P2 * XI1_1 - P2 * XI_11 - P1 * XI_1_1;
}

double Function_T2(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_2(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_2(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_2(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_2(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11 + P2 * XI_11 - P2 * XI1_1 - P1 * XI_1_1;
}

double Function_tau1(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_2(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_2(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_2(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_2(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11 + P2 * XI1_1 - P2 * XI_11 - P1 * XI_1_1;
}

double Function_tau2(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_1(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_1(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_1(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_1(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11 + P2 * XI_11 - P2 * XI1_1 - P1 * XI_1_1;
}

double Function_q1(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_1(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_1(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_1(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_1(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11*XI11 + P2 * XI1_1*XI1_1 + P2 * XI_11*XI_11 + P1 * XI_1_1*XI_1_1;
}

double Function_q2(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_2(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_2(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_2(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_2(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11*XI11 + P2 * XI1_1*XI1_1 + P2 * XI_11*XI_11 + P1 * XI_1_1*XI_1_1;
}

double Function_R(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_12(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_12(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_12(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_12(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11 + P2 * XI1_1 + P2 * XI_11 + P1 * XI_1_1;
}

double Function_r(double z1, double z2, double* param_hat) {
	double XI11 = xiGamma_1(z1, z2, 1., 1., param_hat);
	double XI1_1 = xiGamma_1(z1, z2, 1., -1., param_hat);
	double XI_11 = xiGamma_1(z1, z2, -1., 1., param_hat);
	double XI_1_1 = xiGamma_1(z1, z2, -1., -1., param_hat);
	double _XI11 = xiGamma_2(z1, z2, 1., 1., param_hat);
	double _XI1_1 = xiGamma_2(z1, z2, 1., -1., param_hat);
	double _XI_11 = xiGamma_2(z1, z2, -1., 1., param_hat);
	double _XI_1_1 = xiGamma_2(z1, z2, -1., -1., param_hat);
	double P1 = (1 + q) / 4.;
	double P2 = (1 - q) / 4.;
	return P1 * XI11*_XI11 + P2 * XI1_1*_XI1_1 + P2 * XI_11*_XI_11 + P1 * XI_1_1*_XI_1_1;
}