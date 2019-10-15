#include "OrderParameter.h"

void AllocateNonCavityOrderParam(double** mean1, double** mean2, double** mean12, int Nv) {
	(*mean1) = (double*)calloc(Nv, sizeof(double));
	(*mean2) = (double*)calloc(Nv, sizeof(double));
	(*mean12) = (double*)calloc(Nv, sizeof(double));
}

void CalculateNonCavityOrderParam(double**** u, double**** expt, double* mean1, double* mean2, double* mean12, int Nd, int Nv) {
	for (int i = 0; i < Nv; i++) {
		double expt_i_00 = expt[i][0][0][0] + u[0][i][0][0];
		double expt_i_01 = expt[i][0][0][1] + u[0][i][0][1];
		double expt_i_10 = expt[i][0][1][0] + u[0][i][1][0];
		double expt_i_11 = expt[i][0][1][1] + u[0][i][1][1];
		double denominator = exp(expt_i_10) + exp(expt_i_11) + exp(expt_i_00) + exp(expt_i_01);
		double numerator1 = exp(expt_i_10) + exp(expt_i_11) - exp(expt_i_00) - exp(expt_i_01);
		double numerator2 = exp(expt_i_01) + exp(expt_i_11) - exp(expt_i_00) - exp(expt_i_10);
		double numerator12 = exp(expt_i_00) + exp(expt_i_11) - exp(expt_i_01) - exp(expt_i_10);
		mean1[i] = numerator1 / denominator;
		mean2[i] = numerator2 / denominator;
		mean12[i] = numerator12 / denominator;
	}
}

void CalculateOrderParameter(RBM* rbm, double**** u, double**** expt, double* mean1, double* mean2, double* mean12, double* param) {
	param[T1] = 0.;
	param[T2] = 0.;
	param[tau1] = 0.;
	param[tau2] = 0.;
	param[q1] = 0.;
	param[q2] = 0.;
	param[R] = 0.;
	param[r] = 0.;
	int Nv = rbm->getVisibleNum();
	for (int i = 0; i < Nv; i++) {
		param[T1] += mean1[i] * rbm->weighAt(0, i);
		param[T2]+= mean2[i] * rbm->weighAt(1, i);
		param[tau1] += mean2[i] * rbm->weighAt(0, i);
		param[tau2] += mean1[i] * rbm->weighAt(1, i);
		param[q1] += mean1[i] * mean1[i];
		param[q2] += mean2[i] * mean2[i];
		param[R] += mean12[i];
		param[r] += mean1[i] * mean2[i];
	}
	param[T1] /= Nv;
	param[T2] /= Nv;
	param[tau1] /= Nv;
	param[tau2] /= Nv;
	param[q1] /= Nv;
	param[q2] /= Nv;
	param[R] /= Nv;
	param[r] /= Nv;
}