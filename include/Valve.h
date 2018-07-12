#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>
#include <queue>
#include <functional>
#include "NonPipeComponent.h"
#include "SelfDefinedVariables.h"

using namespace std;

using ValveDataQueue = priority_queue<SignalData, vector<SignalData>, greater<SignalData>>;

class Valve:public NonPipeComponent
{
public:
	Valve();
	Valve(int comp_id, string valve_name, int valve_id, ValveDataQueue valve_status);

public:
	int valve_id;

	vector<int> inlet_connected_comp_set;
	vector<int> outlet_connected_comp_set;
	bool is_reverse;

	double mass_flowrate;
	double inlet_pressure;
	double inlet_temperature;
	double outlet_pressure;
	double outlet_temperature;
	ValveDataQueue status;
};