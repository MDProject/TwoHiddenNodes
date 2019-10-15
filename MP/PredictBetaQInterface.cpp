#include "PredictBetaQInterface.h"

double GradlogZa(RBM* rbm, GaussianParameter* gp, OrderParameter* op, int dataIndex, double** sigma_data) {
	double beta = rbm->beta;
	int Nv = rbm->getVisibleNum();
	double Phi2 = gp->Gamma2[dataIndex][0] + (1 - op->m2[0][dataIndex] * op->m2[0][dataIndex]) / Nv;
	double Phi1 = gp->Gamma1[dataIndex][0] + (1 - op->m1[0][dataIndex] * op->m1[0][dataIndex]) / Nv;
	double Ga1 = gp->G1[dataIndex][0] + sigma_data[dataIndex][0] * op->m1[0][dataIndex] / sqrt(Nv);
	double Ga2 = gp->G2[dataIndex][0] + sigma_data[dataIndex][0] * op->m2[0][dataIndex] / sqrt(Nv);
	double Za = gp->Zeta[dataIndex][0] + (op->q[0][dataIndex] - op->m1[0][dataIndex] * op->m2[0][dataIndex]) / Nv;
	double phi = Za / sqrt(Phi2*Phi1);

	double term1 = beta * Phi2*(1 - phi * phi);
	
	double qs = (gp->Qc[dataIndex][0] + op->q[0][dataIndex] / Nv);
	double term2 = 2 * beta*qs*tanh(beta*beta*qs);

	double sum_phi2 = sqrt(Phi1) + sqrt(Phi2)*phi;
	sum_phi2 = sum_phi2 * sum_phi2;
	double minus_phi2 = sqrt(Phi1) - sqrt(Phi2)*phi;
	minus_phi2 = minus_phi2 * minus_phi2;
	double sum_Ga = Ga1 + Ga2;
	double minus_Ga = Ga1 - Ga2;
	double beta2 = beta * beta;
	double numerator = sum_phi2 * beta*exp(beta2 / 2.*sum_phi2)*cosh(beta*sum_Ga) + exp(beta2 / 2.*sum_phi2)*sum_Ga*sinh(beta*sum_Ga) + beta * minus_phi2*exp(beta2 / 2.*minus_phi2)*cosh(beta*minus_Ga) + exp(beta2 / 2.*minus_phi2)*minus_Ga*sinh(beta*minus_Ga);
	double denominator = exp(beta2 / 2.*sum_phi2)*cosh(beta*(sum_Ga)) + exp(beta2 / 2.*minus_phi2)*cosh(beta*minus_Ga);
	return term1 - term2 + numerator / denominator;
}

void AllocateMessageMatrixGradU(double***** du, int Nd, int Nv) {
	(*du) = (double****)malloc(Nd * sizeof(double***));
	for (int nd = 0; nd < Nd; nd++) {
		(*du)[nd] = (double***)malloc(Nv * sizeof(double**));
		for (int i = 0; i < Nv; i++) {
			(*du)[nd][i] = (double**)malloc(2 * sizeof(double*));
			for (int k = 0; k < 2; k++) {
				(*du)[nd][i][k] = (double*)calloc(2, sizeof(double));
			}
		}
	}
}

void ComputeMessageMatrixGradU(RBM* rbm, double**** du, int Nd, int Nv, double** sigma_data, GaussianParameter* gp) {
	double beta = rbm->beta;
	// 0: -1,-1 [0][0]		1:	-1,1 [0][1]	2:	1,-1 [1][0]	3:	1,1 [1][1], consistent with u
	for (int nd = 0; nd < Nd; nd++) {
		for (int nv = 0; nv < Nv; nv++) {
			double phi = gp->Zeta[nd][nv] / sqrt(gp->Gamma1[nd][nv] * gp->Gamma2[nd][nv]);
			double term1 = beta * gp->Gamma2[nd][nv] * (1. - phi * phi);
			
			double term2_0 = 2.*beta*(gp->Qc[nd][nv] + 1 / Nv)*tanh(beta*beta*(gp->Qc[nd][nv] + 1 / Nv));
			double term2_1 = 2.*beta*(gp->Qc[nd][nv] - 1 / Nv)*tanh(beta*beta*(gp->Qc[nd][nv] - 1 / Nv));
			double term2_2 = term2_1;
			double term2_3 = term2_0;

			double sum_phi2 = sqrt(gp->Gamma1[nd][nv]) + sqrt(gp->Gamma2[nd][nv])*phi;
			sum_phi2 = sum_phi2 * sum_phi2;
			double minus_phi2 = sqrt(gp->Gamma1[nd][nv]) - sqrt(gp->Gamma2[nd][nv])*phi;
			minus_phi2 = minus_phi2 * minus_phi2;
			double inner_sum_term0 = -2.*sigma_data[nd][nv] / sqrt(Nv) + (gp->G1[nd][nv] + gp->G2[nd][nv]);
			double inner_sum_term1 = (gp->G1[nd][nv] + gp->G2[nd][nv]);
			double inner_sum_term2 = inner_sum_term1;
			double inner_sum_term3 = 2.*sigma_data[nd][nv] / sqrt(Nv) + (gp->G1[nd][nv] + gp->G2[nd][nv]);
			double inner_minus_term0 = (gp->G1[nd][nv] - gp->G2[nd][nv]);
			double inner_minus_term3 = inner_minus_term0;
			double inner_minus_term1 = -2.*sigma_data[nd][nv] / sqrt(Nv) + (gp->G1[nd][nv] - gp->G2[nd][nv]);
			double inner_minus_term2 = 2.*sigma_data[nd][nv] / sqrt(Nv) + (gp->G1[nd][nv] - gp->G2[nd][nv]);
			double numerator0 = beta * sum_phi2*exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term0) + exp(beta*beta / 2.*sum_phi2)*inner_sum_term0*sinh(beta*inner_sum_term0) + beta * minus_phi2*exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term0) + exp(beta*beta / 2.*minus_phi2)*inner_minus_term0*sinh(beta*inner_minus_term0);
			double numerator1 = beta * sum_phi2*exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term1) + exp(beta*beta / 2.*sum_phi2)*inner_sum_term1*sinh(beta*inner_sum_term1) + beta * minus_phi2*exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term1) + exp(beta*beta / 2.*minus_phi2)*inner_minus_term1*sinh(beta*inner_minus_term1);
			double numerator2 = beta * sum_phi2*exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term2) + exp(beta*beta / 2.*sum_phi2)*inner_sum_term2*sinh(beta*inner_sum_term2) + beta * minus_phi2*exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term2) + exp(beta*beta / 2.*minus_phi2)*inner_minus_term2*sinh(beta*inner_minus_term2);
			double numerator3 = beta * sum_phi2*exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term3) + exp(beta*beta / 2.*sum_phi2)*inner_sum_term3*sinh(beta*inner_sum_term3) + beta * minus_phi2*exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term3) + exp(beta*beta / 2.*minus_phi2)*inner_minus_term3*sinh(beta*inner_minus_term3);
			double denominator0 = exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term0) + exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term0);
			double denominator1 = exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term1) + exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term1);
			double denominator2 = exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term2) + exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term2);
			double denominator3 = exp(beta*beta / 2.*sum_phi2)*cosh(beta*inner_sum_term3) + exp(beta*beta / 2.*minus_phi2)*cosh(beta*inner_minus_term3);

			du[nd][nv][0][0] = term1 - term2_0 + numerator0 / denominator0;
			du[nd][nv][0][1] = term1 - term2_1 + numerator1 / denominator1;
			du[nd][nv][1][0] = term1 - term2_2 + numerator2 / denominator2;
			du[nd][nv][1][1] = term1 - term2_3 + numerator3 / denominator3;
		}
	}
}

double GradlogZi(RBM* rbm, double**** u, double**** du, double**** expt, int Nd, int vidx) {
	int Nv = rbm->getVisibleNum();
	double dudbeta_0 = 0.;
	double dudbeta_1 = 0.;
	double dudbeta_2 = 0.;
	double dudbeta_3 = 0.;
	// 0: -1,-1 [0][0]		1:	-1,1 [0][1]	2:	1,-1 [1][0]	3:	1,1 [1][1]
	for (int nd = 0; nd < Nd; nd++) {
		dudbeta_0 += du[nd][vidx][0][0];
		dudbeta_1 += du[nd][vidx][0][1];
		dudbeta_2 += du[nd][vidx][1][0];
		dudbeta_3 += du[nd][vidx][1][1];
	}
	// using the shifted exponent tensor 'expt' to build current exponent tensor (make up for the deducted data node)
	double expt0 = expt[vidx][0][0][0] + u[0][vidx][0][0];
	double expt1 = expt[vidx][0][0][1] + u[0][vidx][0][1];
	double expt2 = expt[vidx][0][1][0] + u[0][vidx][1][0];
	double expt3 = expt[vidx][0][1][1] + u[0][vidx][1][1];
	
	double numerator = dudbeta_0 * exp(expt0) + dudbeta_1 * exp(expt1) + dudbeta_2 * exp(expt2) + dudbeta_3 * exp(expt3);
	double denominator = exp(expt0) + exp(expt1) + exp(expt2) + exp(expt3);

	return numerator / denominator;
}

double FeatureCorvariance(MessageParameter* mp, int vidx) {
	double expt0 = mp->expt[vidx][0][0][0] + mp->u[0][vidx][0][0];
	double expt1 = mp->expt[vidx][0][0][1] + mp->u[0][vidx][0][1];
	double expt2 = mp->expt[vidx][0][1][0] + mp->u[0][vidx][1][0];
	double expt3 = mp->expt[vidx][0][1][1] + mp->u[0][vidx][1][1];

	double numerator = exp(expt0) - exp(expt1) - exp(expt2) + exp(expt3);
	double denominator = exp(expt0) + exp(expt1) + exp(expt2) + exp(expt3);

	return numerator / denominator; 
}

