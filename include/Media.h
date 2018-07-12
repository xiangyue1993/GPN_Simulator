#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>

#include "SelfDefinedVariables.h"

using namespace std;

class Media
{
public:
	unordered_map<GasComponentType, double> component;
	enum StateEquationType
	{
		BWRS,
	}equation_type;
	double average_molecular_mass;

public:
	virtual double CalculateDensity(double p, double T) = 0;
	virtual double DpDrho(double rho, double T) = 0;
	virtual double DpDrhoDrho(double rho, double T) = 0;
	virtual double DpDTDT(double rho, double T) = 0;
	virtual double DpDT(double rho, double T) = 0;
	virtual double DrhoDT(double rho, double T) = 0;
	virtual double DpDTDrho(double rho, double T) = 0;
};