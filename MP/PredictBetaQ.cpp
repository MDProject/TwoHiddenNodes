#include "PredictBetaQ.h"

double IterateBeta(RBM* rbm, GaussianParameter* gp, OrderParameter* op, MessageParameter* mp, int Nd, double** sigma_data) {
	int Nv = rbm->getVisibleNum();
	double beta = rbm->beta;
	double Fis = 0.;
	double Fas = 0.;
	for (int i = 0; i < Nv; i++) {
		Fis += GradlogZi(rbm, mp->u, mp->du, mp->expt, Nd, i);
	}
	for (int na = 0; na < Nd; na++) {
		Fas += GradlogZa(rbm, gp, op, na, sigma_data);
	}
	double beta_new = (Fis - (Nv - 1)*Fas) / 2. / Nd;
	return eta * beta + (1 - eta)*beta_new;
}

double IterateQ(MessageParameter* mp, RBM* rbm, double q_odd) {
	int Nv = rbm->getVisibleNum();
	double Q = 0.;
	for (int i = 0; i < Nv; i++) {
		Q += FeatureCorvariance(mp, i);
	}
	Q = Q / (double)Nv;
	return eta * q_odd + (1 - eta)*Q;
}
