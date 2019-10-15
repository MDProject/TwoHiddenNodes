#ifndef DEFINE_ORDERPARAM
#define DEFINE_ORDERPARAM
#include "Parameter.h"
#include "MPInterface.h"

void AllocateNonCavityOrderParam(double** mean1, double** mean2, double** mean12, int Nv);

void CalculateNonCavityOrderParam(double**** u, double**** expt, double* mean1, double* mean2, double* mean12, int Nd, int Nv);

void CalculateOrderParameter(RBM* rbm, double**** u, double**** expt, double* mean1, double* mean2, double* mean12, double* param);

#endif
