#include "RBMInterface.h"
#include "MPInterface.h"
#include "PredictBetaQ.h"
#include "MNIST.h"

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
	std::cout << "Current order param: " << d[i0][j0] << " " << d[i1][j1] << " " << d[i2][j2] << " " << d[i3][j3] << " " << d[i4][j4] << " " << d[i5][j5] << " " << d[i6][j6] << " " << d[i7][j7] << std::endl;
}

double IterationTagPerStep(double** d, int Nd, int Nv) {
	int i0 = 0; int j0 = 0;
	int i1 = floor(Nv / 10); int j1 = floor(Nd / 10);
	int i2 = 2 * floor(Nv / 10); int j2 = 2 * floor(Nd / 10);
	int i3 = 3 * floor(Nv / 10); int j3 = 3 * floor(Nd / 10);
	int i4 = 4 * floor(Nv / 10); int j4 = 4 * floor(Nd / 10);
	int i5 = 5 * floor(Nv / 10); int j5 = 5 * floor(Nd / 10);
	int i6 = 6 * floor(Nv / 10); int j6 = 6 * floor(Nd / 10);
	int i7 = 7 * floor(Nv / 10); int j7 = 7 * floor(Nd / 10);
	return (d[i0][j0] + d[i1][j1] + d[i2][j2] + d[i3][j3] + d[i4][j4] + d[i5][j5] + d[i6][j6] + d[i7][j7]);
}

int main() {
	srand((unsigned)time(NULL));
	int Nv = 28*28;
	int Nh = 2;
	int Nd = 0;
	double q_drive = 0.2; // q true
	RBM* rbm = new RBM(Nh, Nv, 1., true); // beta true 
	// priori of true weight
	GenerateWeightb(rbm, q_drive);
	GenerateHNode(rbm);
	GenerateBiasV(rbm);
	GenerateBiasH(rbm);
	//GibbsSamplingRBM(rbm, 200, 0);

	std::string img_path = "D:\\DeepLearning\\Application\\MNIST_DATA\\MNIST\\raw\\train-images-idx3-ubyte";
	std::string idx_path = "D:\\DeepLearning\\Application\\MNIST_DATA\\MNIST\\raw\\train-labels-idx1-ubyte";
	
	// **********************************
	double** sigma_data;
	double*** img_data;
	unsigned int* label_data;
	unsigned int assign_label[] = { 0,1,2,3 };
	int num_list[] = { 100,100,100,100 };
	Nd = 400;
	
	int Nimg = MnistImageLoader(img_path, &img_data);
	MnistLabelLoader(idx_path, &label_data);
	ExtractAssignedData(assign_label, num_list,2, img_data, label_data, Nimg, 28, 28, &sigma_data);
	CompressData2Binary(sigma_data, Nd, 28 * 28, 125);
	FreeMnistImgMemory(img_data, 60000, 28, 28);
	FreeMnistLabelMemory(label_data);
	// ***********************************
	double** m1, ** m2, ** q;
	InitialOrderParameter(&m1, &m2, &q, Nd, Nv);
	double** Qc, ** G1, **G2, **Gamma1, **Gamma2, **Zeta;
	AllocateCavityQ(Nv, Nd, &Qc);
	AllocateMeanMatrix(Nd, Nv, &G1, &G2);
	AllocateVarianceMatrix(Nd, Nv, &Gamma1, &Gamma2);
	AllocateCovarianceMatrix(Nd, Nv, &Zeta);
	GaussianParameter GP;
	OrderParameter OP;
	OP.m1 = m1;	OP.m2 = m2;	OP.q = q; // initial struct
	MessageParameter MP;
	double**** u, **** expt, **** du;
	AllocateExponentMatrix(Nd, Nv, &expt);
	AllocateMessageMatrixU(rbm, Nd, &u);
	AllocateMessageMatrixGradU(&du, Nd, Nv);
	IterationMonitor(m1, Nd, Nv);

	// A complete iteration loop
	double q_drive_guest = q_drive;
	double beta_guest = 0.8;// rbm->beta;
	double outPrev = -100000;
	for (int n = 0; n < 200; n++) {
		rbm->setBeta(beta_guest);
		double err = 100000.;
		int ac = 0;
		while (abs(err) > epsilon) {
			// compute Gaussian Param
			ComputeCavityQ(q, Nv, Nd, Qc);
			ComputeMeanMatrix(sigma_data, m1, m2, Nd, Nv, G1, G2);
			ComputeVarianceMatrix(m1, m2, Nd, Nv, Gamma1, Gamma2);
			ComputeCovarianceMatrix(m1, m2, q, Nd, Nv, Zeta);
			GP.G1 = G1; GP.G2 = G2; GP.Gamma1 = Gamma1; GP.Gamma2 = Gamma2; GP.Zeta = Zeta; GP.Qc = Qc;

			// compute Message Param
			ComputeMessageMatrixU(rbm, sigma_data, &GP, Nd, u);
			ComputeMessageMatrixGradU(rbm, du, Nd, Nv, sigma_data, &GP);
			ComputeExponentMatrix(Nd, Nv, u, expt, q_drive_guest);
			ShiftExponentMatrixU(Nd, Nv, expt);
			MP.du = du;	MP.expt = expt;	MP.u = u;

			// compute Order Param
			ComputeOrderParameter(m1, m2, q, expt, Nd, Nv);
			OP.m1 = m1;	OP.m2 = m2;	OP.q = q;
			double out = IterationTagPerStep(m1, Nd, Nv);
			err = out - outPrev;
			outPrev = out;
			ac++;
			if (ac > MaxIter) {
				IterationMonitor(m1, Nd, Nv);
				std::cout << "Exceed maximum iteration steps!" << std::endl;
				system("pause");
			}
		}
		IterationMonitor(m1, Nd, Nv);

		beta_guest = IterateBeta(rbm, &GP, &OP, &MP, Nd, sigma_data);
		rbm->setBeta(beta_guest);
		q_drive_guest = IterateQ(&MP, rbm, q_drive_guest);
		std::cout << "Beta: " << beta_guest << "     " << "q:	" << q_drive_guest << std::endl;
	}
	system("pause");
}