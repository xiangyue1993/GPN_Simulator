#pragma once

#include <iostream>
#include <cstdlib>

#include "Media.h"

using namespace std;

#define GRAVITY 9.81
#define R_CONSTANT 8.3143

class BWRS :public Media
{
public:
	BWRS(unordered_map<GasComponentType,double> comp);

public:
	virtual double CalculateDensity(double p, double T);
	virtual double DpDrho(double rho, double T);
	virtual double DpDrhoDrho(double rho, double T);
	virtual double DpDTDT(double rho, double T);
	virtual double DpDT(double rho, double T);
	virtual double DrhoDT(double rho, double T);
	virtual double DpDTDrho(double rho, double T);

private:
	vector<double> BWRS_coefficients;
};