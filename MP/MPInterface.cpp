#include "MPInterface.h"

double Rand() {
	double r = (rand() / (double)RAND_MAX - 0.5)*2.;
	return r;
}

void InitialOrderParameter(double*** m1, double*** m2, double*** q, int Nd, int Nv) {
	(*m1) = (double**)malloc(Nv * sizeof(double*));
	(*m2) = (double**)malloc(Nv * sizeof(double*));
	(*q) = (double**)malloc(Nv * sizeof(double*));
	for (int n = 0; n < Nv; n++) {
		(*m1)[n] = (double*)calloc(Nd, sizeof(double));
		(*m2)[n] = (double*)calloc(Nd, sizeof(double));
		(*q)[n] = (double*)calloc(Nd, sizeof(double));
	}
	srand((unsigned)time(NULL));
	// Initialize the order parameters
	for (int i = 0; i < Nv; i++) {
		for (int j = 0; j < Nd; j++) {
			(*m1)[i][j] = Rand();
			(*m2)[i][j] = Rand();
			(*q)[i][j] = Rand();
			//int a = 1;
		}
	}
}


void AllocateCavityQ(int Nv, int Nd, double*** Qc) {
	(*Qc) = (double**)malloc(Nd * sizeof(double*));
	for (int i = 0; i < Nd; i++) {
		(*Qc)[i] = (double*)calloc(Nv, sizeof(double));
	}
}
void ComputeCavityQ(double** q, int Nv, int Nd, double** Qc) {
	for (int nb = 0; nb < Nd; nb++) {
		double sum_q = 0.;
		for (int j = 0; j < Nv; j++) {
			sum_q += q[j][nb];
		}
		for (int i = 0; i < Nv; i++) {
			Qc[nb][i] = (sum_q - q[i][nb]) / Nv;
		}
	}
}


void AllocateMeanMatrix(int Nd, int Nv, double*** G1, double*** G2) {
	(*G1) = (double**)malloc(Nd * sizeof(double*));
	(*G2) = (double**)malloc(Nd * sizeof(double*));
	for (int n = 0; n < Nd; n++) {
		(*G1)[n] = (double*)calloc(Nv, sizeof(double));
		(*G2)[n] = (double*)calloc(Nv, sizeof(double));
	}
}
void ComputeMeanMatrix(double** sigma_data, double** m1, double** m2, int Nd, int Nv, double** G1, double** G2) {
	for (int nb = 0; nb < Nd; nb++) {
		// for given data b, calculate the summation once
		double sum_sigmab_m1 = 0.;
		double sum_sigmab_m2 = 0.;
		for (int j = 0; j < Nv; j++) {
			sum_sigmab_m1 += sigma_data[nb][j] * m1[j][nb];
			sum_sigmab_m2 += sigma_data[nb][j] * m2[j][nb];
		}
		for (int i = 0; i < Nv; i++) {
			G1[nb][i] = (sum_sigmab_m1 - sigma_data[nb][i] * m1[i][nb]) / sqrt(Nv);
			G2[nb][i] = (sum_sigmab_m2 - sigma_data[nb][i] * m2[i][nb]) / sqrt(Nv);
		}
	}
}

void AllocateVarianceMatrix(int Nd, int Nv, double*** Gamma1, double*** Gamma2) {
	(*Gamma1) = (double**)malloc(Nd * sizeof(double*));
	(*Gamma2) = (double**)malloc(Nd * sizeof(double*));
	for (int n = 0; n < Nd; n++) {
		(*Gamma1)[n] = (double*)calloc(Nv, sizeof(double));
		(*Gamma2)[n] = (double*)calloc(Nv, sizeof(double));
	}
}
void ComputeVarianceMatrix(double** m1, double** m2, int Nd, int Nv, double** Gamma1, double** Gamma2) {
	for (int nb = 0; nb < Nd; nb++) {
		// for given data b, calculate the summation once
		double sum_1_m1 = 0.;
		double sum_1_m2 = 0.;
		for (int j = 0; j < Nv; j++) {
			sum_1_m1 += (1 - m1[j][nb] * m1[j][nb]);
			sum_1_m2 += (1 - m2[j][nb] * m2[j][nb]);
		}
		for (int i = 0; i < Nv; i++) {
			Gamma1[nb][i] = (sum_1_m1 - (1 - m1[i][nb] * m1[i][nb])) / Nv;
			Gamma2[nb][i] = (sum_1_m2 - (1 - m2[i][nb] * m2[i][nb])) / Nv;
		}
	}
}

void AllocateCovarianceMatrix(int Nd, int Nv, double*** Zeta) {
	(*Zeta) = (double**)malloc(Nd * sizeof(double*));
	for (int n = 0; n < Nd; n++) {
		(*Zeta)[n] = (double*)calloc(Nv, sizeof(double));
	}
}
void ComputeCovarianceMatrix(double** m1, double** m2, double**q, int Nd, int Nv, double** Zeta) {
	for (int nb = 0; nb < Nd; nb++) {
		// for given data b, calculate the summation once
		double sum_q_m1m2 = 0.;
		for (int j = 0; j < Nv; j++) {
			sum_q_m1m2 += q[j][nb] - m1[j][nb] * m2[j][nb];
		}
		for (int i = 0; i < Nv; i++) {
			double tmpi = q[i][nb] - m1[i][nb] * m2[i][nb];
			Zeta[nb][i] = (sum_q_m1m2 - tmpi) / Nv;
		}
	}
}

double ComputeMessageU(RBM* rbm, double** sigma_data, GaussianParameter* op, int Nd, int nb, int i,int xi1, int xi2) {
	int Nv = rbm->getVisibleNum();
	double phi = op->Zeta[nb][i] / sqrt(op->Gamma1[nb][i] * op->Gamma2[nb][i]);
	double term1 = rbm->beta*rbm->beta*op->Gamma2[nb][i] * (1 - phi * phi) / 2.;

	double tmp = op->Qc[nb][i] + xi1 * xi2 / Nv;
	double term2 = log(2 * cosh(rbm->beta*rbm->beta*tmp));

	double term3 = rbm->beta*rbm->beta*pow(sqrt(op->Gamma1[nb][i]) + sqrt(op->Gamma2[nb][i])*phi, 2.) / 2.;

	tmp = rbm->beta*op->G1[nb][i] + rbm->beta*op->G2[nb][i] + rbm->beta*sigma_data[nb][i] * (xi1 + xi2) / sqrt(Nv);
	double term4 = log(cosh(tmp));

	double numerator = rbm->beta*(op->G1[nb][i] - op->G2[nb][i] + sigma_data[nb][i] * (xi1 - xi2) / sqrt(Nv));
	double denominator = rbm->beta*(op->G1[nb][i] + op->G2[nb][i] + sigma_data[nb][i] * (xi1 + xi2) / sqrt(Nv));
	double exponent = 0 - 2.*rbm->beta*rbm->beta*op->Zeta[nb][i];
	double term5 = log(1 + exp(exponent)*cosh(numerator) / cosh(denominator));
	
	return term1 - term2 + term3 + term4 + term5;
}

void AllocateMessageMatrixU(RBM* rbm, int Nd, double***** u) {
	int Nv = rbm->getVisibleNum();
	(*u) = (double****)malloc(Nd * sizeof(double***));
	for (int nd = 0; nd < Nd; nd++) {
		(*u)[nd] = (double***)malloc(Nv * sizeof(double**));
		for (int i = 0; i < Nv; i++) {
			(*u)[nd][i] = (double**)malloc(2 * sizeof(double*));
			for (int k = 0; k < 2; k++) {
				(*u)[nd][i][k] = (double*)calloc(2, sizeof(double));
			}
		}
	}
}
void ComputeMessageMatrixU(RBM* rbm, double** sigma_data, GaussianParameter* op, int Nd, double**** u) {
	int Nv = rbm->getVisibleNum();
	for (int nb = 0; nb < Nd; nb++) {
		for (int i = 0; i < Nv; i++) {
			u[nb][i][0][0] = ComputeMessageU(rbm, sigma_data, op, Nd, nb, i, -1, -1);
			u[nb][i][0][1] = ComputeMessageU(rbm, sigma_data, op, Nd, nb, i, -1, 1);
			u[nb][i][1][0] = ComputeMessageU(rbm, sigma_data, op, Nd, nb, i, 1, -1);
			u[nb][i][1][1] = ComputeMessageU(rbm, sigma_data, op, Nd, nb, i, 1, 1);
		}
	}
}

void AllocateExponentMatrix(int Nd, int Nv, double***** expt) {
	(*expt) = (double****)malloc(Nv * sizeof(double***));
	for (int i = 0; i < Nv; i++) {
		(*expt)[i] = (double***)malloc(Nd * sizeof(double**));
		for (int nd = 0; nd < Nd; nd++) {
			(*expt)[i][nd] = (double**)malloc(2 * sizeof(double*));
			for (int k = 0; k < 2; k++) {
				(*expt)[i][nd][k] = (double*)calloc(2, sizeof(double));
			}
		}
	}
}


void ComputeExponentMatrix(int Nd, int Nv, double**** u, double**** expt, double q) {
	for (int i = 0; i < Nv; i++) {
		double sum_u_0_0 = 0.;
		double sum_u_0_1 = 0.;
		double sum_u_1_0 = 0.;
		double sum_u_1_1 = 0.;
		for (int nb = 0; nb < Nd; nb++) {
			sum_u_0_0 += u[nb][i][0][0];
			sum_u_0_1 += u[nb][i][0][1];
			sum_u_1_0 += u[nb][i][1][0];
			sum_u_1_1 += u[nb][i][1][1];
		}
		for (int na = 0; na < Nd; na++) {
			expt[i][na][0][0] = sum_u_0_0 - u[na][i][0][0] + log(priori_probability(-1, -1, q));
			expt[i][na][0][1] = sum_u_0_1 - u[na][i][0][1] + log(priori_probability(-1, 1, q));
			expt[i][na][1][0] = sum_u_1_0 - u[na][i][1][0] + log(priori_probability(1, -1, q));
			expt[i][na][1][1] = sum_u_1_1 - u[na][i][1][1] + log(priori_probability(1, 1, q));
		}
	}
}

void ShiftExponentMatrixU(int Nd, int Nv, double**** expt) {
	for (int i = 0; i < Nv; i++) {
		for (int na = 0; na < Nd; na++) {
			double mean = 0.;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					mean += expt[i][na][k][l];
				}
			}
			mean /= 4;
			expt[i][na][0][0] -= mean;
			expt[i][na][0][1] -= mean;
			expt[i][na][1][0] -= mean;
			expt[i][na][1][1] -= mean;
		}
	}
}

void ComputeOrderParameter(double** m1, double** m2, double** q, double**** expt, int Nd, int Nv) {
	for (int i = 0; i < Nv; i++) {
		for (int na = 0; na < Nd; na++) {
			double numerator = exp(expt[i][na][1][0]) + exp(expt[i][na][1][1]) - exp(expt[i][na][0][1]) - exp(expt[i][na][0][0]);
			double denominator = exp(expt[i][na][1][0]) + exp(expt[i][na][1][1]) + exp(expt[i][na][0][1]) + exp(expt[i][na][0][0]);
			m1[i][na] = numerator / denominator;

			numerator = exp(expt[i][na][0][1]) + exp(expt[i][na][1][1]) - exp(expt[i][na][0][0]) - exp(expt[i][na][1][0]);
			m2[i][na] = numerator / denominator;

			numerator = exp(expt[i][na][1][1]) + exp(expt[i][na][0][0]) - exp(expt[i][na][0][1]) - exp(expt[i][na][1][0]);
			q[i][na] = numerator / denominator;
		}
	}
}



	