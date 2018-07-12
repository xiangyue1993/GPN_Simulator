#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>

#include "SelfDefinedVariables.h"

using namespace std;

class Pipe
{
public:
	Pipe(int id, string name, int inlet_id, int outlet_id, double len, double outerD, 
		double thickness, double roughness, vector<ElevationData> ele, int node_num, 
		vector<InitialPipeData> init_p, vector<InitialPipeData> init_m, vector<InitialPipeData> init_T);

public:
	int pipe_id;
	string pipe_name;
	int inlet_comp_id;
	int outlet_comp_id;

public:
	double length;
	double inter_diameter;
	double outer_diameter;
	double cross_area;
	double wall_thickness;
	double wall_roughness;
	vector<ElevationData> elevation_set;

public:
	double soil_temperature;
	double heat_trans_coe;

public:
	int node_num;
	double dx;

public:
	vector<double> massflowrate;
	vector<double> pressure;
	vector<vector<double>> temperature;
	vector<double> density;
	vector<double> flow_velocity;

public:
	double GetSlope(int node_index);
};
