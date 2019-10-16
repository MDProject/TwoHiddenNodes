#include "RBMInterface.h"
#include "MPInterface.h"
#include "OrderParameter.h"
#include "FreeEnergy.h"

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

void ShowOrderParameter(double* op, int N) {
	for (int i = 0; i < N; i++) {
		std::cout << op[i] << "***";
	}
	std::cout << std::endl;
}

void ShowReferenceParameter(RBM* rbm, double* mean1, double* mean2, double* mean12, int N) {
	std::cout << "xi1" << '\t' << "xi2" << '\t' << "mean 1" << '\t' << "mean 2" << '\t' << "mean12" << std::endl;
	for (int i = 0; i < N; i++) {
		std::cout << rbm->weighAt(0, i) << '\t' << rbm->weighAt(1, i) << '\t' << mean1[i] << '\t' << mean2[i] << '\t' << mean12[i] << std::endl;
	}
}

void ShowFreeEnergy(double f) {
	std::cout << "free energy:\t" << f << std::endl;
}

int main() {
	int Nv = 500;
	int Nh = 2;
	int Nd = Nv*0.25;
	double q_drive = 0.3;
	RBM* rbm = new RBM(Nh, Nv, 1, true);
	// priori of true weight
	GenerateWeightb(rbm, q_drive);
	GenerateHNode(rbm);
	GenerateBiasV(rbm);         
	GenerateBiasH(rbm);
	//GibbsSamplingRBM(rbm, 200, 0);
	
	double** sigma_data;
	//sigma_data = GenerateData(rbm, Nd, 8000, 10, 0);
	GenerateDataSharedMem(rbm, Nd, 10000, 200, 0, &sigma_data, false);
	double** m1, ** m2, ** q;
	InitialOrderParameter(&m1, &m2, &q, Nd, Nv);
	double** Qc, ** G1, **G2, **Gamma1, **Gamma2, **Zeta;
	AllocateCavityQ(Nv, Nd, &Qc);
	AllocateMeanMatrix(Nd, Nv, &G1, &G2);
	AllocateVarianceMatrix(Nd, Nv, &Gamma1, &Gamma2);
	AllocateCovarianceMatrix(Nd, Nv, &Zeta);
	GaussianParameter GP;
	double**** u, **** expt;
	AllocateExponentMatrix(Nd, Nv, &expt);
	AllocateMessageMatrixU(rbm, Nd, &u);
	IterationMonitor(m1, Nd, Nv);
	// A complete iteration loop
	for (int n = 0; n < 70; n++) {
		ComputeCavityQ(q, Nv, Nd, Qc);
		ComputeMeanMatrix(sigma_data, m1, m2, Nd, Nv, G1, G2);
		ComputeVarianceMatrix(m1, m2, Nd, Nv, Gamma1, Gamma2);
		ComputeCovarianceMatrix(m1, m2, q, Nd, Nv, Zeta);

		GP.G1 = G1; GP.G2 = G2; GP.Gamma1 = Gamma1; GP.Gamma2 = Gamma2; GP.Zeta = Zeta; GP.Qc = Qc;

		ComputeMessageMatrixU(rbm, sigma_data, &GP, Nd, u);
		ComputeExponentMatrix(Nd, Nv, u, expt, q_drive);
		ShiftExponentMatrixU(Nd, Nv, expt);
		ComputeOrderParameter(m1, m2, q, expt, Nd, Nv);
		IterationMonitor(m1, Nd, Nv);
	}
	// Compute Order Parameters
	OrderParameter OP;
	MessageParameter MP;
	double Param[8];
	OP.m1 = m1; OP.m2 = m2; OP.q = q;
	MP.u = u; MP.expt = expt;
	double* mean1, *mean2, *mean12;
	AllocateNonCavityOrderParam(&mean1, &mean2, &mean12, Nv);
	CalculateNonCavityOrderParam(u, expt, mean1, mean2, mean12, Nd, Nv);
	CalculateOrderParameter(rbm, u, expt, mean1, mean2, mean12, Param);
	ShowOrderParameter(Param, 8);
	//ShowReferenceParameter(rbm, mean1, mean2, mean12, Nv);
	double fenerg = CalculateFreeEnerg_MP(&MP, &GP, &OP, rbm, sigma_data, q_drive, Nd);
	ShowFreeEnergy(fenerg);
	system("pause");
}  