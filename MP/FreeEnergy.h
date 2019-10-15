#ifndef DEFINE_FREEENERGY
#define DEFINE_FREEENERGY
#include "Priori.h"
#include "MPInterface.h"
#include <math.h>

double CalculateFreeEnerg_MP(MessageParameter* MP, GaussianParameter* GP, OrderParameter* OP, RBM* rbm, double** sigma_data, double q, int Nd);

#endif
