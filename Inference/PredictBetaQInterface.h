#ifndef DEFINE_PREDICT
#define DEFINE_PREDICT

#include "MPInterface.h"
#include "RBMInterface.h"
#include "Parameter.h"
#include <math.h>

void AllocateMessageMatrixGradU(double***** du, int Nd, int Nv);
void ComputeMessageMatrixGradU(RBM* rbm, double**** du, int Nd, int Nv, double** sigma_data, GaussianParameter* gp, OrderParameter* op);

double logZa(RBM* rbm, GaussianParameter* gp, OrderParameter* op, int dataIndex, double** sigma_data);

double logZi(RBM* rbm, double**** u, double q);



#endif
