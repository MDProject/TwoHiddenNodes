#ifndef DEFINE_MPINTERFACE
#define DEFINE_MPINTERFACE

#include <time.h>
#include "RBM.h"
#include "Priori.h"

/*
m1[node index][data index] same for m2 and q
Input parameters are the raw pointers 
*/

class GaussianParameter {
public:
	double** G1, ** G2, ** Gamma1, ** Gamma2, ** Zeta, ** Qc;
};

class OrderParameter {
public:
	double** m1, **m2, **q;
};

void InitialOrderParameter(double*** m1, double*** m2, double*** q, int Nd, int Nv);

// Qc[data index][visible nodes index]
void AllocateCavityQ(int Nv, int Nd, double*** Qc);
void ComputeCavityQ(double** q, int Nv, int Nd, double** Qc);

// Nd: data number	Nv: visible nodes number	G1 and G2 are 2D matrix of cavity mean 
void AllocateMeanMatrix(int Nd, int Nv, double*** G1, double*** G2);
void ComputeMeanMatrix(double** sigma_data, double** m1, double** m2, int Nd, int Nv, double** G1, double** G2);

void AllocateVarianceMatrix(int Nd, int Nv, double*** Gamma1, double*** Gamma2);
void ComputeVarianceMatrix(double** m1, double** m2, int Nd, int Nv, double** Gamma1, double** Gamma2);

void AllocateCovarianceMatrix(int Nd, int Nv, double*** Zeta);
void ComputeCovarianceMatrix(double** m1, double** m2, double**q, int Nd, int Nv, double** Zeta);

double ComputeMessageU(RBM* rbm, double** sigma_data, GaussianParameter* op, int Nd, int nb, int i, int xi1, int xi2);

// Here u is a 4th-order tensor pointer [data index][visible nodes index][+][-]
void AllocateMessageMatrixU(RBM* rbm, int Nd, double***** u);
void ComputeMessageMatrixU(RBM* rbm, double** sigma_data, GaussianParameter* op, int Nd, double**** u);

// Exponent expt is a 4-th order tensor pointer [visible nodes index][data index][+][-]
void AllocateExponentMatrix(int Nd, int Nv, double***** expt);
void ComputeExponentMatrix(int Nd, int Nv, double**** u, double**** expt, double q);

// Shifting the exponents to avoid the exploring or vanishing exponential values
void ShiftExponentMatrixU(int Nd, int Nv, double**** expt);

void ComputeOrderParameter(double** m1, double** m2, double** q, double**** expt, int Nd, int Nv);

#endif
