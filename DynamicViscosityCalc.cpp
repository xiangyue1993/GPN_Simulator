#include "DynamicViscosityCalc.h"

double DynamicViscosityCalc(double rho, double T) {
	double re_rho = rho / 1.206;

	double c1 = 2.57 + 0.2781 * re_rho + 1063.6 / T;
	double c2 = 1.11 + 0.04 * c1;
	double c3 = 2.415*(7.77 + 0.1844*re_rho)*pow(T, 1.5)*1e-4 / (122.4 + 377.58*re_rho + 1.8*T);

	return c3 * exp(c1*pow(rho / 1000, c2)) * 1e-3;
}