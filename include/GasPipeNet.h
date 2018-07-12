#pragma once

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unordered_map>
#include <map>

#include "NonPipeComponent.h"
#include "Valve.h"
#include "Compressor.h"
#include "GasSource.h"
#include "Pipe.h"
#include "Media.h"
#include "BWRS.h"

#include "json\json.h"

using namespace std;

class GasPipeNet
{
public:
	GasPipeNet();
	GasPipeNet(Json::Value& jsonFileContent);
public:
	map<pair<int, int>, int> net_structure; //Do not use unordered_map for pair type key, unless you have an effective hash solution for it
	vector<NonPipeComponent*> non_pipe_component_set;
	vector<Pipe*> pipe_set;
	Media* media;

	int pipe_num;
	int valve_num;
	int compressor_num;
	int gas_source_num;
	int connection_node_num;

public:
	void Output();
};