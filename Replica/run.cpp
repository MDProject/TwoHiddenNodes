#include "Function.h"
#include "Integrator.h"
#include <iostream>
#include "Define.h"
#include "FreeEnergy.h"

// recursion start from [T1,T2,tau1,tau2,q1,q2,R,r] 

void PrintInfo(double* param, double* param_hat) {
	for (int i = 0; i < 8; i++) {
		std::cout << param[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < 8; i++) {
		std::cout << param_hat[i] << " ";
	}
	std::cout << "---------------------------------" << std::endl;
}

int InitialParam(double* param) {
	bool ifLoop = true;
	int ac = 0;
	while (ifLoop) {
		for (int i = 0; i < 8; i++) {
			param[i] = Rand();
		}
		if (B(param) >= 0 && (C(param) >= 0)) {
			ifLoop = false;
		}
		ac++;
	}
	return ac;
}

void ShowFreeEnergy(double f) {
	std::cout << "free energy:\t" << f << std::endl;
}

int main() {
	double Param[8] = { 0.697737, 0.698987, 0.667203, 0.669043, 0.70419, 0.703036, 0.539728, 0.685023 };
	double Param_hat[8] = { 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 };
	int N = 20;
	int loop = 300000;
	srand((unsigned)time(NULL));
	int ns = InitialParam(Param);
	std::cout << "Loop " << ns << "times to get initial order parameters" << std::endl;
	for (int i = 0; i < N; i++) {
		Param_hat[0] = INT4Gaussian(Function_T1_hat, Param, loop);
		
		Param_hat[1] = INT4Gaussian(Function_T2_hat, Param, loop);
		Param_hat[2] = INT4Gaussian(Function_tau1_hat, Param, loop);
		Param_hat[3] = INT4Gaussian(Function_tau2_hat, Param, loop);
		Param_hat[4] = INT4Gaussian(Function_q1_hat, Param, loop);
		Param_hat[5] = INT4Gaussian(Function_q2_hat, Param, loop);
		Param_hat[6] = INT4Gaussian(Function_R_hat, Param, loop) - alpha * beta*beta*tanh(beta*beta*Param[R]);
		Param_hat[7] = INT4Gaussian(Function_r_hat, Param, loop);

		//Param_hat[tau1] = 0.;
		//Param_hat[tau2] = 0.;
		/*Param_hat[T1] = 0.;
		Param_hat[T2] = 0.;
		Param_hat[R] = 0.;
		Param_hat[r] = 0.;*/

		Param[0] = INT2Gaussian(Function_T1, Param_hat, loop);
		if (isnan(Param[0])) {
			std::cout << std::endl;
		}
		Param[1] = INT2Gaussian(Function_T2, Param_hat, loop);
		Param[2] = INT2Gaussian(Function_tau1, Param_hat, loop);
		Param[3] = INT2Gaussian(Function_tau2, Param_hat, loop);
		Param[4] = INT2Gaussian(Function_q1, Param_hat, loop);
		Param[5] = INT2Gaussian(Function_q2, Param_hat, loop);
		Param[6] = INT2Gaussian(Function_R, Param_hat, loop);
		Param[7] = INT2Gaussian(Function_r, Param_hat, loop);

		/*Param[T1] = 0.;
		Param[T2] = 0.;
		Param[R] = 0.;
		Param[r] = 0.;*/
		PrintInfo(Param, Param_hat);
	}
	double free_energy = CalculateFreeEnergy(Param, Param_hat, loop);
	ShowFreeEnergy(free_energy);
	system("pause");
}