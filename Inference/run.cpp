#include "RBMInterface.h"
#include "MPInterface.h"

void IterationMonitor(double** d, int Nd, int Nv) {
	// observe the behavior of order parameter 
	int i0 = 0; int j0 = 0;
	int i1 = floor(Nv / 10); int j1 = floor(Nd / 10);
	int i2 = 2 * floor(Nv / 10); int j2 = 2 * floor(Nd / 10);
	int i3 = 3 * floor(Nv / 10); int j3 = 3 * floor(Nd / 10);
	int i4 = 4 * floor(Nv / 10); int j4 = 4 * floor(Nd / 10);
	int i5 = 5 * floor(Nv / 10); int j5 = 5 * floor(Nd / 10);
	int i6 = 6 * floor(Nv / 10); int j6 = 6 * floor(Nd / 10);
	int i7 = 7 * floor(Nv / 10); int j7 = 7 * floor(Nd / 10);
	std::cout << d[i0][j0] << " " << d[i1][j1] << " " << d[i2][j2] << " " << d[i3][j3] << " " << d[i4][j4] << " " << d[i5][j5] << " " << d[i6][j6] << " " << d[i7][j7] << std::endl;
}

int main() {
	int Nv = 900;
	int Nh = 2;
	int Nd = 200;
	double q_drive = 0.2;
	RBM* rbm = new RBM(Nh, Nv, 1, true);
	// priori of true weight
	GenerateWeightb(rbm, q_drive);
	GenerateHNode(rbm);
	GenerateBiasV(rbm);         
	GenerateBiasH(rbm);
	//GibbsSamplingRBM(rbm, 200, 0);
	
	double** sigma_data = GenerateData(rbm, Nd, 8000, 10, 0);
	double** m1, ** m2, ** q;
	InitialOrderParameter(&m1, &m2, &q, Nd, Nv);
	double** Qc, ** G1, **G2, **Gamma1, **Gamma2, **Zeta;
	AllocateCavityQ(Nv, Nd, &Qc);
	AllocateMeanMatrix(Nd, Nv, &G1, &G2);
	AllocateVarianceMatrix(Nd, Nv, &Gamma1, &Gamma2);
	AllocateCovarianceMatrix(Nd, Nv, &Zeta);
	GaussianParameter OP;
	double**** u, **** expt;
	AllocateExponentMatrix(Nd, Nv, &expt);
	AllocateMessageMatrixU(rbm, Nd, &u);
	IterationMonitor(m1, Nd, Nv);
	// A complete iteration loop
	for (int n = 0; n < 200; n++) {
		ComputeCavityQ(q, Nv, Nd, Qc);
		ComputeMeanMatrix(sigma_data, m1, m2, Nd, Nv, G1, G2);
		ComputeVarianceMatrix(m1, m2, Nd, Nv, Gamma1, Gamma2);
		ComputeCovarianceMatrix(m1, m2, q, Nd, Nv, Zeta);

		OP.G1 = G1; OP.G2 = G2; OP.Gamma1 = Gamma1; OP.Gamma2 = Gamma2; OP.Zeta = Zeta; OP.Qc = Qc;

		ComputeMessageMatrixU(rbm, sigma_data, &OP, Nd, u);
		ComputeExponentMatrix(Nd, Nv, u, expt, q_drive);
		ShiftExponentMatrixU(Nd, Nv, expt);
		ComputeOrderParameter(m1, m2, q, expt, Nd, Nv);
		IterationMonitor(m1, Nd, Nv);
	}
	system("pause");
}  