#include "GibbsSampling.h"

void GibbsSamplingRBM(RBM* rbm, int NStep, int tag) {
	int Nv = rbm->getVisibleNum();
	int Nh = rbm->getHiddenNum();
	int N = NStep;
	if (tag == 1) {// start from visible layer (data)
		for (int k = 0; k < Nh; k++) {
			double tmp = 0;
			for (int j = 0; j < Nv; j++) {
				tmp += rbm->weighAt(k, j)*rbm->vnode[j];
			}
			tmp /= sqrt(Nv);
			tmp += rbm->b[k];
			double numerator = exp(0 - tmp);
			double denominator = 2 * cosh(tmp);
			double Prob = numerator / denominator;
			double r = rand() / (double)RAND_MAX;
			if (r < Prob) {
				rbm->hnode[k] = 1;
			}
			else {
				rbm->hnode[k] = -1;
			}
		}
		N -= 1;
	}
	for (int n = 0; n < N; n++) {
		// hidden -- visible
		for (int k = 0; k < Nv; k++) {
			double tmp = 0.;
			for (int i = 0; i < Nh; i++) {
				tmp += rbm->weighAt(i, k)*rbm->hnode[i];
			}
			tmp /= sqrt(Nv);
			tmp += rbm->a[k];
			double numerator = exp(0 - tmp);
			double denominator = 2 * cosh(tmp);
			double Prob = numerator / denominator;
			double r = rand() / (double)RAND_MAX;
			if (r < Prob) {
				rbm->vnode[k] = 1;
			}
			else {
				rbm->vnode[k] = -1;
			}
		}
		// visible -- hidden
		for (int k = 0; k < Nh; k++) {
			double tmp = 0;
			for (int j = 0; j < Nv; j++) {
				tmp += rbm->weighAt(k, j)*rbm->vnode[j];
			}
			tmp /= sqrt(Nv);
			tmp += rbm->b[k];
			double numerator = exp(0 - tmp);
			double denominator = 2 * cosh(tmp);
			double Prob = numerator / denominator;
			double r = rand() / (double)RAND_MAX;
			if (r < Prob) {
				rbm->hnode[k] = 1;
			}
			else {
				rbm->hnode[k] = -1;
			}
		}
	}
}
