#ifndef DEFINE_PREDICTINTERFACE
#define DEFINE_PREDICTINTERFACE

#include "MPInterface.h"
#include "RBMInterface.h"
#include "Parameter.h"
#include <math.h>

void AllocateMessageMatrixGradU(double***** du, int Nd, int Nv);
void ComputeMessageMatrixGradU(RBM* rbm, double**** du, int Nd, int Nv, double** sigma_data, GaussianParameter* gp);

double GradlogZa(RBM* rbm, GaussianParameter* gp, OrderParameter* op, int dataIndex, double** sigma_data);

double GradlogZi(RBM* rbm, double**** u, double**** du, double**** expt, int Nd, int vidx);

double FeatureCorvariance(MessageParameter* mp, int vidx);

#endif

