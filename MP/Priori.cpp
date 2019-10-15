#include "Priori.h"
#include <math.h>

double priori_probability(int xi1, int xi2, double q) {
	if (xi1*xi2 == 1) {
		return (1 + q) / 4.;
	}
	else {
		return (1 - q) / 4.;
	}
}

double priori_probability_(int xi1, int xi2, double q) {
	double J0 = atan(q);
	return (double)exp(J0*xi1*xi2) / 4. / cosh(J0);
}