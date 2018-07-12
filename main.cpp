#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "GasPipeNet.h"
#include "Solve.h"

#include "json\json.h"

using namespace std;

int main() {
	Json::CharReaderBuilder builder;
	builder["collectComments"] = false;
	JSONCPP_STRING errs;
	Json::Value root;

	ifstream jsonFile("Input.json");
	if (Json::parseFromStream(builder, jsonFile, &root, &errs)) {
		double duration = root["SimulationConfig"]["duration"].asDouble();
		double time_increment = root["SimulationConfig"]["time_increment"].asDouble();

		GasPipeNet* _net = new GasPipeNet(root);
		Solve(_net, duration, time_increment);
	}
	
	return 0;
}