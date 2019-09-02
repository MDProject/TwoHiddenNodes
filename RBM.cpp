#include "RBM.h"

RBM::RBM(int N1, int N2, double beta, bool ifBinary) :Nh(N1), N(N2), beta(beta){
	// indicate weight property
	this->ifBinary = ifBinary;
	AllocateMem();
}

RBM::~RBM() {
	free(a);
	free(b);
	free(vnode);
	free(hnode);
	for (int i = 0; i < Nh; i++) {
		if (ifBinary) {
			free(weightb[i]);
		}
		else {
			free(weight[i]);
		}
	}
	if (ifBinary) {
		free(weightb);
	}
	else {
		free(weight);
	}
}

void RBM::AllocateMem() {
	// allocate weight memory
	if (ifBinary) {
		weight = NULL;
		weightb = (int**)malloc(sizeof(int*)*Nh);
		for (int i = 0; i < Nh; i++) {
			weightb[i] = (int*)malloc(sizeof(int)*N);
		}
	}
	else {
		weight = (double**)malloc(sizeof(double*)*Nh); 
		for (int i = 0; i < Nh; i++) {
			weight[i] = (double*)malloc(sizeof(double)*N);
		}
		weightb = NULL;
	}
	// allocate bias memory
	a = (double*)malloc(sizeof(double)*N);
	vnode = (double*)malloc(sizeof(double)*N);
	b = (double*)malloc(sizeof(double)*Nh);
	hnode = (double*)malloc(sizeof(double)*Nh);
}

double RBM::weighAt(int i, int j) {
	if (ifBinary) {
		return weightb[i][j];
	}
	else {
		return weight[i][j];
	}
}

void RBM::saveWeightb(int i, int j, int w) {
	weightb[i][j] = w;
}

void RBM::saveWeight(int i, int j, double w) {
	weight[i][j] = w;
}