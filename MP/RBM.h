#ifndef DEFINE_RBM
#define DEFINE_RBM

#include "Parameter.h"

class RBM {
	int Nh;
	int N;
	// binary weight is independently realized here
	double** weight;
	int** weightb; // weight[hidden node index][visible node index]
	bool ifBinary;
	void AllocateMem();
public:
	double* a; // visible bias
	double* b; // hidden bias
	double* vnode; // only store one copy of data
	double* hnode;
	double beta;
	RBM(int N1, int N2, double beta, bool ifBinary = false);
	~RBM();
	void setBeta(double b);
	int getVisibleNum() { return N; }
	int getHiddenNum() { return Nh; }
	double weighAt(int i, int j);
	void saveWeightb(int i, int j, int w);
	void saveWeight(int i, int j, double w);
};


#endif
