#pragma once

#include <cstdlib>
#include <iostream>

#include "NonPipeComponent.h"

using namespace std;

class Compressor :public NonPipeComponent
{
public:
	int comp_id;

	vector<int> inlet_connected_comp_set;
	vector<int> outlet_connected_comp_set;
	bool is_reverse;

	vector<double> mass_flowrate;
	vector<double> inlet_pressure;
	vector<double> inlet_temperature;
	vector<double> outlet_pressure;
	vector<double> outlet_temperature;
	double pressure_ratio;
};