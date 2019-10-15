#include "Priori.h"

double probability_priori(double xi1, double xi2, double q) {
	if (xi1*xi2 > 0) {
		return (1 + q) / 4.;
	}
	else {
		return (1 - q) / 4.;
	}
}