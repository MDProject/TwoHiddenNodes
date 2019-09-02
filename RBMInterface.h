#ifndef DEFINE_RBMINTERFACE
#define DEFINE_RBMINTERFACE

#include "RBM.h"
#include "GibbsSampling.h"

void GenerateWeightb(RBM* rbm, double q);
void GenerateHNode(RBM* rbm);
void GenerateVNode(RBM* rbm);
void GenerateBiasV(RBM* rbm);
void GenerateBiasH(RBM* rbm);

// data stored in visible nodes of RBM	[data index][visible nodes index]
double** GenerateData(RBM* rbm, int Nd, int NStep, int outFreq, int tag);

#endif

