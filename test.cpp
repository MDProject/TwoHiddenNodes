
#include "RBMInterface.h"
#include "MPInterface.h"

int main() {
	int Nv = 50;
	int Nh = 2;
	int Nd = 20;
	RBM* rbm = new RBM(Nh, Nv, 1, true);
	// priori of true weight
	GenerateWeightb(rbm, 0.2);
	GenerateHNode(rbm);
	GenerateBiasV(rbm);         
	GenerateBiasH(rbm);
	//GibbsSamplingRBM(rbm, 200, 0);
	
	double** sigma_data = GenerateData(rbm, Nd, 200, 10, 0);
	double** m1, ** m2, ** q;
	InitialOrderParameter(&m1, &m2, &q, Nd, Nv);
	
	// A complete iteration loop

	double** Qc, ** G1, **G2, **Gamma1, **Gamma2, **Zeta;
	ComputeCavityQ(q, Nv, Nd, &Qc);
	ComputeMeanMatrix(sigma_data, m1, m2, Nd, Nv, &G1, &G2);
	ComputeVarianceMatrix(m1, m2, Nd, Nv, &Gamma1, &Gamma2);
	ComputeCovarianceMatrix(m1, m2, q, Nd, Nv, &Zeta);

	OrderParameter OP;
	OP.G1 = G1; OP.G2 = G2; OP.Gamma1 = Gamma1; OP.Gamma2 = Gamma2; OP.Zeta = Zeta; OP.Qc = Qc;
	double**** u, **** expt;
	ComputeMessageMatrixU(rbm, sigma_data, &OP, Nd, &u);
	ComputeExponentMatrix(Nd, Nv, u, &expt);
	ShiftExponentMatrixU(Nd, Nv, expt);
	ComputeOrderParameter(m1, m2, q, expt, Nd, Nv);

	system("pause");
}  