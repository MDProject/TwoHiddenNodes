#ifndef FUNCTION_H
#define FUNCTION_H


double B(double* param);
double C(double* param);


double Function_T1_hat(double t, double x, double u, double u_prime, double* param);
double Function_T2_hat(double t, double x, double u, double u_prime, double* param);
double Function_tau1_hat(double t, double x, double u, double u_prime, double* param);
double Function_tau2_hat(double t, double x, double u, double u_prime, double* param);
double Function_q1_hat(double t, double x, double u, double u_prime, double* param);
double Function_q2_hat(double t, double x, double u, double u_prime, double* param);
double Function_R_hat(double t, double x, double u, double u_prime, double* param);
double Function_r_hat(double t, double x, double u, double u_prime, double* param);

double Function_T1(double z1, double z2, double* param_hat);
double Function_T2(double z1, double z2, double* param_hat);
double Function_tau1(double z1, double z2, double* param_hat);
double Function_tau2(double z1, double z2, double* param_hat);
double Function_q1(double z1, double z2, double* param_hat);
double Function_q2(double z1, double z2, double* param_hat);
double Function_R(double z1, double z2, double* param_hat);
double Function_r(double z1, double z2, double* param_hat);


#endif
