#include "RBMInterface.h"

int rand1() {
	double r = rand() / (double)RAND_MAX;
	if (r > 0.5) {
		return 1;
	}
	else {
		return -1;
	}
}

void GenerateWeightb(RBM* rbm, double q) {
	int Nv = rbm->getVisibleNum();
	//int Nh = rbm->getHiddenNum();
	for (int j = 0; j < Nv; j++) {
		rbm->saveWeightb(0, j, rand1());
		double ratio = 0.5*(1 + q);
		if (rand() / (double)RAND_MAX < ratio) {
			rbm->saveWeightb(1, j, rbm->weighAt(0, j));
		}
		else {
			rbm->saveWeightb(1, j, -1 / rbm->weighAt(0, j));
			int y = 0;
		}
	}
}

void GenerateHNode(RBM* rbm) {
	int Nh = rbm->getHiddenNum();
	for (int i = 0; i < Nh; i++) {
		rbm->hnode[i] = rand1();
	}
}

void GenerateVNode(RBM* rbm) {
	int Nv = rbm->getVisibleNum();
	for (int i = 0; i < Nv; i++) {
		rbm->vnode[i] = rand1();
	}
}

void GenerateBiasV(RBM* rbm) {
	int Nv = rbm->getVisibleNum();
	for (int i = 0; i < Nv; i++) {
		rbm->a[i] = 0;
	}
}

void GenerateBiasH(RBM* rbm) {
	int Nh = rbm->getHiddenNum();
	for (int i = 0; i < Nh; i++) {
		rbm->b[i] = 0;
	}
}

double** GenerateData(RBM* rbm, int Nd, int NStep, int outFreq, int tag) {
	int Nv = rbm->getVisibleNum();
	double** dataArray = (double**)calloc(Nd, sizeof(double*));
	for (int n = 0; n < Nd; n++) {
		dataArray[n] = (double*)calloc(Nv, sizeof(double));
	}
	GibbsSamplingRBM(rbm, NStep, tag);
	for (int n = 0; n < Nd; n++) {
		GibbsSamplingRBM(rbm, outFreq, tag);
		memcpy(dataArray[n], rbm->vnode, Nv * sizeof(double));
	}
	return dataArray;
}

void GenerateDataSharedMem(RBM* rbm, int Nd, int NStep, int outFreq, int tag, double*** sigma_data, bool ifValid) {
	int Nv = rbm->getVisibleNum();
	if (!ifValid) {
		(*sigma_data) = (double**)calloc(Nd, sizeof(double*));
		for (int n = 0; n < Nd; n++) {
			(*sigma_data)[n] = (double*)calloc(Nv, sizeof(double));
		}
	}
	GibbsSamplingRBM(rbm, NStep, tag);
	for (int n = 0; n < Nd; n++) {
		GibbsSamplingRBM(rbm, outFreq, tag);
		memcpy((*sigma_data)[n], rbm->vnode, Nv * sizeof(double));
	}
}