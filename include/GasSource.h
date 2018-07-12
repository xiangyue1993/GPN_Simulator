#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>
#include <queue>
#include <functional>

#include "NonPipeComponent.h"
#include "SelfDefinedVariables.h"

using namespace std;

using SourceDataQueue = priority_queue<SignalData, vector<SignalData>, greater<SignalData>>;

class GasSource:public NonPipeComponent
{
public:
	GasSource(int comp_id, string source_name, int source_id, BoundaryConditionType bound_type, SourceDataQueue& p, SourceDataQueue& m, SourceDataQueue& T);

public:
	int source_id;
	BoundaryConditionType bound_type;

	vector<int> connected_comp_set;

	SourceDataQueue pressure;
	SourceDataQueue mass_flowrate;
	SourceDataQueue temperature;
	double density;
};