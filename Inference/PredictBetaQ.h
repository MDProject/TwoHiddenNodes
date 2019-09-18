#ifndef DEFINE_PREDICT
#define DEFINE_PREDICT

#include "PredictBetaQInterface.h"

double IterateBeta(RBM* rbm, GaussianParameter* gp, OrderParameter* op, MessageParameter* mp, int Nd, double** sigma_data);

double IterateQ(MessageParameter* mp, RBM* rbm);

#endif
