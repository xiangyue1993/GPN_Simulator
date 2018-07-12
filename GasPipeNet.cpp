#include "GasPipeNet.h"

static int output_index = 0;

unordered_map<string, GasComponentType> nameToindex = {
	{"CH4", GasComponentType::CH4},
	{"C2H6", GasComponentType::C2H6 },
	{ "C3H8", GasComponentType::C3H8 },
	{ "i_C4H10", GasComponentType::i_C4H10 },
	{ "n_C4H10", GasComponentType::n_C4H10 },
	{ "i_C5H12", GasComponentType::i_C5H12 },
	{ "n_C5H12", GasComponentType::n_C5H12 },
	{ "C6H14", GasComponentType::C6H14 },
	{ "C7H16", GasComponentType::C7H16 },
	{ "C8H18", GasComponentType::C8H18 },
	{ "C9H20", GasComponentType::C9H20 },
	{ "C10H22", GasComponentType::C10H22 },
	{ "C2H4", GasComponentType::C2H4 },
	{ "C3H6", GasComponentType::C3H6 },
	{ "N2", GasComponentType::N2 },
	{ "CO2", GasComponentType::CO2 },
	{ "H2S", GasComponentType::H2S },
	{ "H2O", GasComponentType::H2O },
};

bool CmpFuncforNonpipeSet(const NonPipeComponent* a, const NonPipeComponent* b) {
	return a->nonpipe_comp_id < b->nonpipe_comp_id;
}

bool CmpFuncForPipeSet(const Pipe* a, const Pipe* b) {
	return a->pipe_id < b->pipe_id;
}

GasPipeNet::GasPipeNet() {
	
}

GasPipeNet::GasPipeNet(Json::Value& root) {
	pipe_num = 0;
	valve_num = 0;
	compressor_num = 0;
	gas_source_num = 0;
	pipe_set.clear();
	non_pipe_component_set.clear();
	media = NULL;

	unordered_map<GasComponentType, double> gas_comp;
	auto members = root["GasComponent"].getMemberNames();
	for (auto iter = members.begin(); iter != members.end(); iter++) {
		GasComponentType cur_comptype = nameToindex[*iter];
		gas_comp[cur_comptype] = root["GasComponent"][*iter].asDouble() / 100;
	}
	if (root["SimulationConfig"]["StateEquationType"].asString() == "BWRS") {
		media = new BWRS(gas_comp);
	}

	for (auto iter1 = root["NonPipeComponent"].begin(); iter1 != root["NonPipeComponent"].end(); iter1++) {
		Json::Value elements = root["NonPipeComponent"][iter1.name()];

		if (elements["class_name"].asString() == "Gas_Source") {
			int comp_id = elements["comp_id"].asInt();
			int source_id = elements["source_id"].asInt();

			bool isPressureSpecified= false;
			bool isMassflowrateSpecified = false;
			SourceDataQueue p, m, T;
			BoundaryConditionType boundtype = None;
			for (auto iter2 = elements["boundary_conditions"].begin(); iter2 != elements["boundary_conditions"].end(); iter2++) {
				string mem_name = iter2.name();
				for (int i = 0; i < elements["boundary_conditions"][mem_name].size(); i++) {
					string signal_mode = elements["boundary_conditions"][mem_name][i][0].asString();
					double time_occur = elements["boundary_conditions"][mem_name][i][1].asDouble();
					double data_value = elements["boundary_conditions"][mem_name][i][2].asDouble();
						
					SignalData temp(signal_mode, time_occur, data_value);

					if (mem_name == "pressure") {
						isPressureSpecified = true;
						p.push(temp);
					}
						
					if (mem_name == "massflowrate") {
						isMassflowrateSpecified = true;
						m.push(temp);
					}

					if (mem_name == "temperature") {
						T.push(temp);
					}
				}
			}
			if(isPressureSpecified && isMassflowrateSpecified){
				boundtype = BothPressureAndMassFlowRate;
			}
			else if (isPressureSpecified)
			{
				boundtype = OnlyPressure;
			}
			else if (isMassflowrateSpecified)
			{
				boundtype = OnlyMassFlowRate;
			}

			NonPipeComponent* p_source = new GasSource(comp_id, iter1.name(), source_id, boundtype, p, m, T);
			non_pipe_component_set.push_back(p_source);
			gas_source_num++;
		}

		if (elements["class_name"].asString() == "Valve") {
			int comp_id = elements["comp_id"].asInt();
			int valve_id = elements["valve_id"].asInt();

			ValveDataQueue status;
			for (int i = 0; i < elements["status"].size(); i++) {
				string status_mode = elements["status"][i][0].asString();
				double cur_time = elements["status"][i][1].asDouble();
				if (status_mode == "on") {
					double pressure_ratio = elements["status"][i][2].asDouble();
					status.push(SignalData(status_mode, cur_time, pressure_ratio));
				}
				else {
					status.push(SignalData(status_mode, cur_time, 0));
				}
			}

			NonPipeComponent* p_valve = new Valve(comp_id, iter1.name(), valve_id, status);
			non_pipe_component_set.push_back(p_valve);
			valve_num++;
		}
	}
	sort(non_pipe_component_set.begin(), non_pipe_component_set.end(), CmpFuncforNonpipeSet);

	for (auto iter1 = root["PipeComponent"].begin(); iter1 != root["PipeComponent"].end(); iter1++) {
		Json::Value elements = root["PipeComponent"][iter1.name()];
		int pipe_id = elements["pipe_id"].asInt();
		double length = elements["length"].asDouble();
		double outer_diameter = elements["outer_diameter"].asDouble();
		double wall_thickness = elements["wall_thickness"].asDouble();
		double wall_roughness = elements["wall_roughness"].asDouble();

		vector<ElevationData> ele;
		for (int i = 0; i < elements["elevation"].size(); i++) {
			double horizontal_distance = elements["elevation"][i][0].asDouble();
			double elevation = elements["elevation"][i][1].asDouble();
			ele.push_back(ElevationData(horizontal_distance, elevation));
		}

		int node_num = elements["node_num"].asInt();
		vector<InitialPipeData> init_p, init_m, init_T;
		for (auto iter2 = elements["initial_conditions"].begin(); iter2 != elements["initial_conditions"].end(); iter2++) {
			string mem_name = iter2.name();
			for (int i = 0; i < iter2->size(); i++) {
				double horizontal_distance = elements["initial_conditions"][mem_name][i][0].asDouble();
				double data_value = elements["initial_conditions"][mem_name][i][1].asDouble();

				if (mem_name == "pressure") {
					init_p.push_back(InitialPipeData(horizontal_distance, data_value));
				}
				if (mem_name == "massflowrate") {
					init_m.push_back(InitialPipeData(horizontal_distance, data_value));
				}
				if (mem_name == "temperature") {
					init_T.push_back(InitialPipeData(horizontal_distance, data_value));
				}
			}
		}

		int inlet_comp_id = elements["inlet_comp_id"].asInt();
		int outlet_comp_id = elements["outlet_comp_id"].asInt();
		net_structure[make_pair(inlet_comp_id, outlet_comp_id)] = pipe_id;

		if (non_pipe_component_set[inlet_comp_id]->comp_type == GasSource_Type) {
			GasSource* p_pre = (GasSource *)non_pipe_component_set[inlet_comp_id];
			p_pre->connected_comp_set.push_back(outlet_comp_id);
		}
		else if (non_pipe_component_set[inlet_comp_id]->comp_type == Valve_Type) {
			Valve* p_pre = (Valve *)non_pipe_component_set[inlet_comp_id];
			p_pre->outlet_connected_comp_set.push_back(outlet_comp_id);
		}
		else if (non_pipe_component_set[inlet_comp_id]->comp_type == Compressor_Type) {
			Compressor* p_pre = (Compressor *)non_pipe_component_set[inlet_comp_id];
			p_pre->outlet_connected_comp_set.push_back(outlet_comp_id);
		}
		
		if (non_pipe_component_set[outlet_comp_id]->comp_type == GasSource_Type) {
			GasSource* p_next = (GasSource *)non_pipe_component_set[outlet_comp_id];
			p_next->connected_comp_set.push_back(inlet_comp_id);
		}
		else if(non_pipe_component_set[outlet_comp_id]->comp_type == Valve_Type)
		{
			Valve* p_next = (Valve *)non_pipe_component_set[outlet_comp_id];
			p_next->inlet_connected_comp_set.push_back(inlet_comp_id);
		}
		else if (non_pipe_component_set[outlet_comp_id]->comp_type == Compressor_Type) {
			Compressor* p_next = (Compressor *)non_pipe_component_set[outlet_comp_id];
			p_next->inlet_connected_comp_set.push_back(inlet_comp_id);
		}

		Pipe* p_pipe = new Pipe(pipe_id, iter1.name(), inlet_comp_id, outlet_comp_id, length, outer_diameter,
			wall_thickness, wall_roughness, ele, node_num, init_p, init_m, init_T);
		for (int i = 0; i < p_pipe->node_num; i++) {
			double cur_pressure = p_pipe->pressure[i];
			double cur_temperature = p_pipe->temperature[1][i];
			double cur_massflowrate = p_pipe->massflowrate[i];
			p_pipe->density.push_back(media->CalculateDensity(cur_pressure, cur_temperature));
			p_pipe->flow_velocity.push_back(cur_massflowrate / p_pipe->cross_area / p_pipe->density[i]);
		}
		pipe_set.push_back(p_pipe);
		pipe_num++;
	}

	sort(pipe_set.begin(), pipe_set.end(), CmpFuncForPipeSet);
}

void GasPipeNet::Output() {
	Json::Value net_data;
	
	double cur_time = 0;
	for (int i = 0; i < non_pipe_component_set.size(); i++) {
		if (non_pipe_component_set[i]->comp_type == GasSource_Type) {
			GasSource* p_source = (GasSource *)non_pipe_component_set[i];
			cur_time = p_source->pressure.top().time_occur;

			Json::Value source_data;
			source_data["pressure"] = Json::Value(p_source->pressure.top().value);
			source_data["massflowrate"] = Json::Value(p_source->mass_flowrate.top().value);
			source_data["density"] = Json::Value(p_source->density);
			source_data["temperature"] = Json::Value(p_source->temperature.top().value);

			string source_name = p_source->comp_name;
			net_data["NonPipeComponent"][source_name] = source_data;
		}

		if (non_pipe_component_set[i]->comp_type == Valve_Type) {
			Valve* p_valve = (Valve *)non_pipe_component_set[i];

			Json::Value valve_data;
			valve_data["massflowrate"] = Json::Value(p_valve->mass_flowrate);
			valve_data["inlet_pressure"] = Json::Value(p_valve->inlet_pressure);
			valve_data["outlet_pressure"] = Json::Value(p_valve->outlet_pressure);
			valve_data["inlet_temperature"] = Json::Value(p_valve->inlet_temperature);
			valve_data["outlet_temperature"] = Json::Value(p_valve->outlet_temperature);

			string valve_name = p_valve->comp_name;
			net_data["NonPipeComponent"][valve_name] = valve_data;
		}
	}

	for (int i = 0; i < pipe_set.size(); i++) {
		Pipe* p_pipe = pipe_set[i];

		Json::Value pipe_data;
		for (int j = 0; j < p_pipe->node_num; j++) {
			Json::Value item_p, item_m, item_T, item_v, item_d;
			double pos = p_pipe->dx * j;

			item_p.append(Json::Value(pos));
			item_p.append(Json::Value(p_pipe->pressure[j]));
			pipe_data["pressure"].append(item_p);
			item_m.append(Json::Value(pos));
			item_m.append(Json::Value(p_pipe->massflowrate[j]));
			pipe_data["massflowrate"].append(item_m);
			item_v.append(Json::Value(pos));
			item_v.append(Json::Value(p_pipe->flow_velocity[j]));
			pipe_data["velocity"].append(item_v);
			item_d.append(Json::Value(pos));
			item_d.append(Json::Value(p_pipe->density[j]));
			pipe_data["density"].append(item_d);
			item_T.append(Json::Value(pos));
			item_T.append(Json::Value(p_pipe->temperature[1][j]));
			pipe_data["temperature"].append(item_T);
		}

		string pipe_name = p_pipe->pipe_name;
		net_data["PipeComponent"][pipe_name] = pipe_data;
	}

	net_data["current_time"] = Json::Value(cur_time);


	string file_name = "output/Output_" + to_string(output_index) + ".json";

	ofstream ofs;
	ofs.open(file_name);
	ofs << net_data.toStyledString();
	ofs.close();
	output_index++;

	return;
}