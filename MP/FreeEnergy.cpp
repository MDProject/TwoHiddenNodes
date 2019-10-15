#include "FreeEnergy.h"

// Substract the mean from exponents
double CalculateLnZi(double**** u, double q, int Nd, int Nv, int ni) {
	double sum_u_0_0 = 0.;
	double sum_u_0_1 = 0.;
	double sum_u_1_0 = 0.;
	double sum_u_1_1 = 0.;
	for (int nb = 0; nb < Nd; nb++) {
		sum_u_0_0 += u[nb][ni][0][0];
		sum_u_0_1 += u[nb][ni][0][1];
		sum_u_1_0 += u[nb][ni][1][0];
		sum_u_1_1 += u[nb][ni][1][1];
	}
	// add priori
	sum_u_0_0 += log(priori_probability(-1, -1, q));
	sum_u_0_1 += log(priori_probability(-1, 1, q));
	sum_u_1_0 += log(priori_probability(1, -1, q));
	sum_u_1_1 += log(priori_probability(1, 1, q));
	double mu = (sum_u_0_0 + sum_u_0_1 + sum_u_1_0 + sum_u_1_1) / 4.;
	sum_u_0_0 -= mu;
	sum_u_0_1 -= mu;
	sum_u_1_0 -= mu;
	sum_u_1_1 -= mu;
	double lnZi = log(exp(sum_u_0_0) + exp(sum_u_0_1) + exp(sum_u_1_0) + exp(sum_u_1_1)) + mu;
	return lnZi;
}

double CalculateLnZa(GaussianParameter* GP, OrderParameter* OP, RBM* rbm, double** sigma_data, int Nd, int na) {
	int Nv = rbm->getVisibleNum();
	double beta = rbm->beta;
	// Non Cavity variables
	double zeta = GP->Zeta[na][0] + (OP->q[0][na] - OP->m1[0][na] * OP->m2[0][na]) / Nv;
	double gamma1 = GP->Gamma1[na][0] + (1. - OP->m1[0][na] * OP->m1[0][na]) / Nv;
	double gamma2 = GP->Gamma2[na][0] + (1. - OP->m2[0][na] * OP->m2[0][na]) / Nv;
	double phi = zeta / sqrt(gamma1*gamma2);
	double Qa = GP->Qc[na][0] + OP->q[0][na] / Nv;
	double Ga1 = GP->G1[na][0] + sigma_data[na][0] * OP->m1[0][na] / sqrt(Nv);
	double Ga2 = GP->G2[na][0] + sigma_data[na][0] * OP->m2[0][na] / sqrt(Nv);

	double term1 = rbm->beta*rbm->beta*gamma2*(1. - phi * phi) / 2.;
	double term2 = log(2.*cosh(beta*beta*Qa));
	double term3 = beta * beta*(sqrt(gamma1) + sqrt(gamma2)*phi)*(sqrt(gamma1) + sqrt(gamma2)*phi) / 2.;
	double term4 = log(cosh(beta*(Ga1 + Ga2)));
	double numerator = cosh(beta*(Ga1 - Ga2));
	double denominator = cosh(beta*(Ga1 + Ga2));
	double term5 = log(1 + exp(0.-2.*beta*beta*zeta)*numerator / denominator);
	double lnZa = term1 - term2 + term3 + term4 + term5;
	return lnZa;
}

double CalculateFreeEnerg_MP(MessageParameter* MP, GaussianParameter* GP, OrderParameter* OP, RBM* rbm, double** sigma_data, double q, int Nd) {
	double fi = 0.;
	int Nv = rbm->getVisibleNum();
	for (int i = 0; i < Nv; i++) {
		fi += CalculateLnZi(MP->u, q, Nd, Nv, i);
	}
	double fa = 0.;
	for (int na = 0; na < Nd; na++) {
		fa += CalculateLnZa(GP, OP, rbm, sigma_data, Nd, na);
	}
	double f = (fi - (Nv - 1)*fa) / Nv;
	return f;
}