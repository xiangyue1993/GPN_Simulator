#include "Solve.h"
#include "DynamicViscosityCalc.h"
#include "FrictionCoefficientCalc.h"

#include <armadillo>
#include <iomanip>

int time_index;
double dt;
unordered_map<int, arma::mat> fundamental_solution_sets;

void UpdateBound(GasPipeNet* _net);
void Hydraulic_Solve(GasPipeNet* _net);
void GeneratingMatrixAndVector(GasPipeNet* _net, arma::sp_mat& connect_equation_coe_matrix, arma::vec& connect_equation_cons_vec);
void BackSubstitution(GasPipeNet* _net, arma::vec& connect_equation_unknown_vec);
arma::mat Solve_Pipe(GasPipeNet* _net, Pipe* cur_pipe);

void Solve(GasPipeNet* _net, double duration, double time_increment) {
	time_index = 0;
	dt = time_increment;
	fundamental_solution_sets.clear();

	while (dt * time_index < duration) {
		time_index++;
		UpdateBound(_net);
		Hydraulic_Solve(_net);
		_net->Output();
		cout << "Current time: " << setiosflags(ios::fixed) <<setprecision(2) << dt * time_index << " s" << endl;
	}
}

void UpdateBound(GasPipeNet* _net) {
	for (int i = 0; i < _net->non_pipe_component_set.size(); i++) {
		if (_net->non_pipe_component_set[i]->comp_type == GasSource_Type) {
			GasSource* p_source = (GasSource *)(_net->non_pipe_component_set[i]);
			if (p_source->bound_type == OnlyPressure || p_source->bound_type == BothPressureAndMassFlowRate) {
				SignalData temp(p_source->pressure.top());
				while (!p_source->pressure.empty() && time_index * dt >= p_source->pressure.top().time_occur) { //"While" loop here for considering case that the time between two adjacent signal might be less than time increment, i.e. time_ocurr[i+1] - time_ocurr[i] > dt
					temp = p_source->pressure.top();
					p_source->pressure.pop();
				}
				if (!p_source->pressure.empty() && p_source->pressure.top().signal_mode == "linear") {
					double value_dif = p_source->pressure.top().value - temp.value;
					double time_dif = p_source->pressure.top().time_occur - temp.time_occur;
					double change_rate = value_dif / time_dif;

					double cur_value = temp.value + change_rate * (time_index * dt - temp.time_occur);
					p_source->pressure.push(SignalData("linear", time_index * dt, cur_value));
				}
				else {
					double cur_value = temp.value;
					p_source->pressure.push(SignalData("linear", time_index * dt, cur_value));
				}
			}

			if (p_source->bound_type == OnlyMassFlowRate || p_source->bound_type == BothPressureAndMassFlowRate) {
				SignalData temp(p_source->mass_flowrate.top());
				while (!p_source->mass_flowrate.empty() && time_index * dt >= p_source->mass_flowrate.top().time_occur) {
					temp = p_source->mass_flowrate.top();
					p_source->mass_flowrate.pop();
				}
				if (!p_source->mass_flowrate.empty() && p_source->mass_flowrate.top().signal_mode == "linear") {
					double value_dif = p_source->mass_flowrate.top().value - temp.value;
					double time_dif = p_source->mass_flowrate.top().time_occur - temp.time_occur;
					double change_rate = value_dif / time_dif;

					double cur_value = temp.value + change_rate * (time_index * dt - temp.time_occur);
					p_source->mass_flowrate.push(SignalData("linear", time_index * dt, cur_value));
				}
				else {
					double cur_value = temp.value;
					p_source->mass_flowrate.push(SignalData("linear", time_index * dt, cur_value));
				}
			}

			SignalData temp(p_source->temperature.top());
			while (!p_source->temperature.empty() && time_index * dt >= p_source->temperature.top().time_occur) {
				temp = p_source->temperature.top();
				p_source->temperature.pop();
			}
			if (!p_source->temperature.empty() && p_source->temperature.top().signal_mode == "linear") {
				double value_dif = p_source->temperature.top().value - temp.value;
				double time_dif = p_source->temperature.top().time_occur - temp.time_occur;
				double change_rate = value_dif / time_dif;

				double cur_value = temp.value + change_rate * (time_index * dt - temp.time_occur);
				p_source->temperature.push(SignalData("linear", time_index * dt, cur_value));
			}
			else {
				double cur_value = temp.value;
				p_source->temperature.push(SignalData("linear", time_index * dt, cur_value));
			}
		}

		if (_net->non_pipe_component_set[i]->comp_type == Valve_Type) {
			Valve* p_valve = (Valve *)_net->non_pipe_component_set[i];

			SignalData temp(p_valve->status.top());
			while (!p_valve->status.empty() && time_index * dt >= p_valve->status.top().time_occur) {
				temp = p_valve->status.top();
				p_valve->status.pop();
			}
			temp.time_occur = time_index * dt;
			p_valve->status.push(temp);
		}
	}
}

void Hydraulic_Solve(GasPipeNet* _net) {
	int matrix_order = 4 * _net->pipe_num + 4 * _net->compressor_num + 4 * _net->valve_num + 2 * _net->gas_source_num;
	arma::sp_mat connect_equation_coe_matrix(matrix_order, matrix_order);
	arma::vec connect_equation_cons_vec(matrix_order, arma::fill::zeros);

	GeneratingMatrixAndVector(_net, connect_equation_coe_matrix, connect_equation_cons_vec);

	//arma::mat temp(connect_equation_coe_matrix);
	//temp.print("coe_matrix");
	//connect_equation_cons_vec.print("cons_vec");
	arma::vec connect_equation_unknown_vec = arma::spsolve(connect_equation_coe_matrix, connect_equation_cons_vec);
	//connect_equation_unknown_vec.print("unknown");

	BackSubstitution(_net, connect_equation_unknown_vec);

	return;
}

void GeneratingMatrixAndVector(GasPipeNet* _net, arma::sp_mat& connect_equation_coe_matrix, arma::vec& connect_equation_cons_vec) {
	int equation_index = 0;
	fundamental_solution_sets.clear();

	for (int i = 0; i < _net->pipe_num; i++) {
		arma::mat fundamental_solution = Solve_Pipe(_net, _net->pipe_set[i]);
		//fundamental_solution.print("fund_solution");
		int row_num = fundamental_solution.n_rows;
		int pipe_index = _net->pipe_set[i]->pipe_id;
		int start_unknown_index = 4 * pipe_index;

		connect_equation_coe_matrix(equation_index, start_unknown_index) = -fundamental_solution(1, 0);
		connect_equation_coe_matrix(equation_index, start_unknown_index + 1) = 1;
		connect_equation_coe_matrix(equation_index, start_unknown_index + 3) = -fundamental_solution(1, 1);
		connect_equation_cons_vec(equation_index) = fundamental_solution(1, 2);
		equation_index++;

		connect_equation_coe_matrix(equation_index, start_unknown_index) = -fundamental_solution(row_num - 2, 0);
		connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
		connect_equation_coe_matrix(equation_index, start_unknown_index + 3) = -fundamental_solution(row_num - 2, 1);
		connect_equation_cons_vec(equation_index) = fundamental_solution(row_num - 2, 2);
		equation_index++;

		fundamental_solution_sets[pipe_index] = fundamental_solution;
	}

	for (int i = 0; i < _net->non_pipe_component_set.size(); i++) {
		if (_net->non_pipe_component_set[i]->comp_type == Compressor_Type) {
			Compressor* p_compressor = (Compressor *)(_net->non_pipe_component_set[i]);
			int comp_index = p_compressor->comp_id;
			int start_unknown_index = 4 * _net->pipe_num + 4 * comp_index;

			if (!p_compressor->is_reverse) {
				connect_equation_coe_matrix(equation_index, start_unknown_index) = -p_compressor->pressure_ratio;
				connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
			}
			else {
				connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
				connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = -p_compressor->pressure_ratio;
			}
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			connect_equation_coe_matrix(equation_index, start_unknown_index + 1) = 1;
			connect_equation_coe_matrix(equation_index, start_unknown_index + 3) = -1;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			arma::rowvec temp1(connect_equation_cons_vec.size(), arma::fill::zeros);
			temp1(start_unknown_index + 1) = -1;
			for (int j = 0; j < p_compressor->inlet_connected_comp_set.size(); j++) {//inlet
				int noncomp1_id = p_compressor->nonpipe_comp_id;
				int noncomp2_id = p_compressor->inlet_connected_comp_set[j];
				if (_net->net_structure.find(make_pair(noncomp1_id, noncomp2_id)) != _net->net_structure.end()) {//flow out
					int pipe_index = _net->net_structure.find(make_pair(noncomp1_id, noncomp2_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp1(pipe_index * 4 + 1) = -1;
				}
				else {//flow in
					int pipe_index = _net->net_structure.find(make_pair(noncomp2_id, noncomp1_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4 + 2) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp1(pipe_index * 4 + 3) = 1;
				}
			}
			connect_equation_coe_matrix.row(equation_index) = temp1;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			arma::rowvec temp2(connect_equation_cons_vec.size(), arma::fill::zeros);
			temp2(start_unknown_index + 3) = 1;
			for (int j = 0; j < p_compressor->outlet_connected_comp_set.size(); j++) {//outlet
				int noncomp1_id = p_compressor->nonpipe_comp_id;
				int noncomp2_id = p_compressor->outlet_connected_comp_set[j];
				if (_net->net_structure.find(make_pair(noncomp1_id, noncomp2_id)) != _net->net_structure.end()) {//flow out
					int pipe_index = _net->net_structure.find(make_pair(noncomp1_id, noncomp2_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp2(pipe_index * 4 + 1) = -1;
				}
				else {//flow in
					int pipe_index = _net->net_structure.find(make_pair(noncomp2_id, noncomp1_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4 + 2) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp2(pipe_index * 4 + 3) = 1;
				}
			}
			connect_equation_coe_matrix.row(equation_index) = temp2;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;
		}

		if (_net->non_pipe_component_set[i]->comp_type == Valve_Type) {
			Valve* p_valve = (Valve*)(_net->non_pipe_component_set[i]);
			int valve_index = p_valve->valve_id;
			int start_unknown_index = 4 * _net->pipe_num + 4 * _net->compressor_num + 4 * valve_index;

			if (p_valve->status.top().signal_mode == "on") {
				if (!p_valve->is_reverse) {
					connect_equation_coe_matrix(equation_index, start_unknown_index) = -p_valve->status.top().value;
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
				}
				else {
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = -p_valve->status.top().value;
				}
			}
			else {
				connect_equation_coe_matrix(equation_index, start_unknown_index + 1) = 1;
			}
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			connect_equation_coe_matrix(equation_index, start_unknown_index + 1) = 1;
			connect_equation_coe_matrix(equation_index, start_unknown_index + 3) = -1;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			arma::rowvec temp1(connect_equation_cons_vec.size(), arma::fill::zeros);
			temp1(start_unknown_index + 1) = -1;
			for (int j = 0; j < p_valve->inlet_connected_comp_set.size(); j++) {//inlet
				int noncomp1_id = p_valve->nonpipe_comp_id;
				int noncomp2_id = p_valve->inlet_connected_comp_set[j];
				if (_net->net_structure.find(make_pair(noncomp1_id, noncomp2_id)) != _net->net_structure.end()) {//flow out
					int pipe_index = _net->net_structure.find(make_pair(noncomp1_id, noncomp2_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp1(pipe_index * 4 + 1) = -1;
				}
				else {//flow in
					int pipe_index = _net->net_structure.find(make_pair(noncomp2_id, noncomp1_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4 + 2) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp1(pipe_index * 4 + 3) = 1;
				}
			}
			connect_equation_coe_matrix.row(equation_index) = temp1;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;

			arma::rowvec temp2(connect_equation_cons_vec.size(), arma::fill::zeros);
			temp2(start_unknown_index + 3) = 1;
			for (int j = 0; j < p_valve->outlet_connected_comp_set.size(); j++) {//outlet
				int noncomp1_id = p_valve->nonpipe_comp_id;
				int noncomp2_id = p_valve->outlet_connected_comp_set[j];
				if (_net->net_structure.find(make_pair(noncomp1_id, noncomp2_id)) != _net->net_structure.end()) {//flow out
					int pipe_index = _net->net_structure.find(make_pair(noncomp1_id, noncomp2_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp2(pipe_index * 4 + 1) = -1;
				}
				else {//flow in
					int pipe_index = _net->net_structure.find(make_pair(noncomp2_id, noncomp1_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index + 2) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4 + 2) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp2(pipe_index * 4 + 3) = 1;
				}
			}
			connect_equation_coe_matrix.row(equation_index) = temp2;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;
		}

		if (_net->non_pipe_component_set[i]->comp_type == GasSource_Type) {
			GasSource* p_source = (GasSource*)(_net->non_pipe_component_set[i]);
			int source_index = p_source->source_id;
			int start_unknown_index = 4 * _net->pipe_num + 4 * _net->compressor_num + 4 * _net->valve_num + 2 * source_index;

			if (p_source->bound_type == OnlyPressure || p_source->bound_type == BothPressureAndMassFlowRate) {
				connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
				connect_equation_cons_vec(equation_index) = p_source->pressure.top().value;
				equation_index++;
			}

			if (p_source->bound_type == OnlyMassFlowRate || p_source->bound_type == BothPressureAndMassFlowRate) {
				connect_equation_coe_matrix(equation_index, start_unknown_index + 1) = 1;
				connect_equation_cons_vec(equation_index) = p_source->mass_flowrate.top().value;
				equation_index++;
			}

			arma::rowvec temp(connect_equation_cons_vec.size(), arma::fill::zeros);
			temp(start_unknown_index + 1) = 1;
			for (int j = 0; j < p_source->connected_comp_set.size(); j++) {
				int noncomp1_id = p_source->nonpipe_comp_id;
				int noncomp2_id = p_source->connected_comp_set[j];
				if (_net->net_structure.find(make_pair(noncomp1_id, noncomp2_id)) != _net->net_structure.end()) {//flow out
					int pipe_index = _net->net_structure.find(make_pair(noncomp1_id, noncomp2_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp(pipe_index * 4 + 1) = -1;
				}
				else {//flow in
					int pipe_index = _net->net_structure.find(make_pair(noncomp2_id, noncomp1_id))->second;
					Pipe* p_curpipe = _net->pipe_set[pipe_index];
					connect_equation_coe_matrix(equation_index, start_unknown_index) = 1;
					connect_equation_coe_matrix(equation_index, pipe_index * 4 + 2) = -1;
					connect_equation_cons_vec(equation_index) = 0;
					equation_index++;
					temp(pipe_index * 4 + 3) = 1;
				}
			}
			connect_equation_coe_matrix.row(equation_index) = temp;
			connect_equation_cons_vec(equation_index) = 0;
			equation_index++;
		}
	}

	return;
}

void BackSubstitution(GasPipeNet* _net, arma::vec& connect_equation_unknown_vec) {
	for (int i = 0; i < _net->pipe_num; i++) {
		Pipe* p_curpipe = _net->pipe_set[i];
		int start_unknown_index = 4 * p_curpipe->pipe_id;
		vector<double> new_pressure(p_curpipe->node_num);
		vector<double> new_massflowrate(p_curpipe->node_num);
		vector<double> new_density(p_curpipe->node_num);
		vector<double> new_velocity(p_curpipe->node_num);

		arma::mat fundamental_solution = fundamental_solution_sets[p_curpipe->pipe_id];
		//fundamental_solution.print("fund_solution");
		double p1 = connect_equation_unknown_vec(start_unknown_index);
		double mn = connect_equation_unknown_vec(start_unknown_index + 3);
		for (int j = 0; j < p_curpipe->node_num; j++) {
			new_pressure[j] = fundamental_solution(j * 2, 0) * p1 + fundamental_solution(j * 2, 1) * mn + fundamental_solution(j * 2, 2);
			new_massflowrate[j] = fundamental_solution(j * 2 + 1, 0)*p1 + fundamental_solution(j * 2 + 1, 1)*mn + fundamental_solution(j * 2 + 1, 2);
			new_density[j] = _net->media->CalculateDensity(new_pressure[j], p_curpipe->temperature[1][j]);
			new_velocity[j] = new_massflowrate[j] / (new_density[j] * p_curpipe->cross_area);
		}

		p_curpipe->pressure = new_pressure;
		p_curpipe->massflowrate= new_massflowrate;
		p_curpipe->density= new_density;
		p_curpipe->flow_velocity = new_velocity;
		p_curpipe->temperature.erase(p_curpipe->temperature.begin());
		p_curpipe->temperature.push_back(p_curpipe->temperature[0]);
	}

	for (int i = 0; i < _net->non_pipe_component_set.size(); i++) {
		if (_net->non_pipe_component_set[i]->comp_type == Compressor_Type) {
			Compressor* p_compressor = (Compressor*)(_net->non_pipe_component_set[i]);
			int start_unknown_index = 4 * _net->pipe_num + 4 * p_compressor->comp_id;

			p_compressor->inlet_pressure.push_back(connect_equation_unknown_vec(start_unknown_index));
			p_compressor->mass_flowrate.push_back(connect_equation_unknown_vec(start_unknown_index + 1));
			p_compressor->outlet_pressure.push_back(connect_equation_unknown_vec(start_unknown_index + 2));
		}

		if (_net->non_pipe_component_set[i]->comp_type == Valve_Type) {
			Valve* p_valve = (Valve*)(_net->non_pipe_component_set[i]);
			int start_unknown_index = 4 * _net->pipe_num + 4 * _net->compressor_num + 4 * p_valve->valve_id;

			p_valve->inlet_pressure = connect_equation_unknown_vec(start_unknown_index);
			p_valve->mass_flowrate = connect_equation_unknown_vec(start_unknown_index + 1);
			p_valve->outlet_pressure = connect_equation_unknown_vec(start_unknown_index + 2);

			if (p_valve->mass_flowrate >= 0) {
				p_valve->is_reverse = false;
			}
			else {
				p_valve->is_reverse = true;
			}
		}

		if (_net->non_pipe_component_set[i]->comp_type == GasSource_Type) {
			GasSource* p_source = (GasSource*)(_net->non_pipe_component_set[i]);
			int start_unknown_index = 4 * _net->pipe_num + 4 * _net->compressor_num + 4 * _net->valve_num + 2 * p_source->source_id;

			if (p_source->bound_type == OnlyMassFlowRate || p_source->bound_type == None) {
				if (!p_source->pressure.empty()) p_source->pressure.pop();

				double new_value = connect_equation_unknown_vec(start_unknown_index);
				p_source->pressure.push(SignalData("linear", time_index * dt, new_value));
			}

			if (p_source->bound_type == OnlyPressure || p_source->bound_type == None) {
				if (!p_source->mass_flowrate.empty()) p_source->mass_flowrate.pop();

				double new_value = connect_equation_unknown_vec(start_unknown_index + 1);
				p_source->mass_flowrate.push(SignalData("linear", time_index * dt, new_value));
			}

			double p = p_source->pressure.top().value;
			double T = p_source->temperature.top().value;
			p_source->density = _net->media->CalculateDensity(p, T);
		}
	}

	return;
}

arma::mat Solve_Pipe(GasPipeNet* _net, Pipe* cur_pipe) {
	int matrix_order = cur_pipe->node_num * 2;
	arma::sp_mat pipe_equation_coe_matrix(matrix_order, matrix_order);
	arma::mat pipe_equation_cons_matrix(matrix_order, 3, arma::fill::zeros);

	pipe_equation_coe_matrix(0, 0) = 1;
	pipe_equation_cons_matrix(0, 0) = 1;
	pipe_equation_coe_matrix(matrix_order - 1, matrix_order - 1) = 1;
	pipe_equation_cons_matrix(matrix_order - 1, 1) = 1;

	for (int i = 1; i < matrix_order - 1; i++) {
		int cur_node_index = (i - 1) / 2;
		double mid_p = (cur_pipe->pressure[cur_node_index + 1] + cur_pipe->pressure[cur_node_index]) / 2;
		double mid_m = (cur_pipe->massflowrate[cur_node_index + 1] + cur_pipe->massflowrate[cur_node_index]) / 2;
		double mid_T = (cur_pipe->temperature[1][cur_node_index + 1] + cur_pipe->temperature[1][cur_node_index]) / 2;
		double mid_old_T = (cur_pipe->temperature[0][cur_node_index + 1] + cur_pipe->temperature[0][cur_node_index]) / 2;
		//double mid_rho = (cur_pipe->density[cur_node_index + 1] + cur_pipe->density[cur_node_index]) / 2;
		double mid_rho = _net->media->CalculateDensity(mid_p, mid_T);
		double mid_v = mid_m / mid_rho / cur_pipe->cross_area;
		double slope = cur_pipe->GetSlope(cur_node_index);
		double cross_area = cur_pipe->cross_area;
		double dx = cur_pipe->dx;
		double diameter = cur_pipe->inter_diameter;

		if (i % 2) {
			double dp_drho = _net->media->DpDrho(mid_rho, mid_T);
			double ddp_drho_drho = _net->media->DpDrhoDrho(mid_rho, mid_T);
			double dm_dx = (cur_pipe->massflowrate[cur_node_index + 1] - cur_pipe->massflowrate[cur_node_index]) / dx;
			double dT_dt = (mid_T - mid_old_T) / dt;
			double dp_dT = _net->media->DpDT(mid_rho, mid_T);
			double ddp_dT_drho = _net->media->DpDTDrho(mid_rho, mid_T);

			double B11 = 0;
			double B12 = 1 / cross_area * dp_drho;
			double F1 = dp_dT * dT_dt;
			double G11 = 1 / cross_area * ddp_drho_drho / dp_drho * dm_dx;
			double G12 = 0;
			double dF_dU_11 = dT_dt * ddp_dT_drho / dp_drho;
			double dF_dU_12 = 0;

			double CE11 = 0.5 / dt - 1 / dx * B11 + 0.5 * (G11 - dF_dU_11);
			double CE12 = 0 / dt - 1 / dx * B12 + 0.5 * (G12 - dF_dU_12);
			double DW11 = 0.5 / dt + 1 / dx * B11 + 0.5 * (G11 - dF_dU_11);
			double DW12 = 0 / dt + 1 / dx * B12 + 0.5 * (G12 - dF_dU_12);
			double H1 = F1 + (G11 - dF_dU_11 + 1 / dt) * mid_p + (G12 - dF_dU_12 + 0)*mid_m;

			pipe_equation_coe_matrix(i, cur_node_index * 2) = CE11;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 1) = CE12;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 2) = DW11;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 3) = DW12;
			pipe_equation_cons_matrix(i, 0) = 0;
			pipe_equation_cons_matrix(i, 1) = 0;
			pipe_equation_cons_matrix(i, 2) = H1;
		}
		else {
			double dp_drho = _net->media->DpDrho(mid_rho, mid_T);
			double drho_dp = 1 / dp_drho;
			double dy_vis = DynamicViscosityCalc(mid_rho, mid_T);
			double friction_coe = FrictionCoefficientCalc(cur_pipe, mid_rho, mid_v, dy_vis);
			double dm_dx = (cur_pipe->massflowrate[cur_node_index + 1] - cur_pipe->massflowrate[cur_node_index]) / dx;
			double dp_dx = (cur_pipe->pressure[cur_node_index + 1] - cur_pipe->pressure[cur_node_index]) / dx;
			double dT_dx = (cur_pipe->temperature[1][cur_node_index + 1] - cur_pipe->temperature[1][cur_node_index]) / dx;
			double drho_dT = _net->media->DrhoDT(mid_rho, mid_T);
			double ddp_drho_dp = _net->media->DpDrhoDrho(mid_rho, mid_T) / dp_drho;
			double ddrho_dp_dp = -ddp_drho_dp / (dp_drho*dp_drho);
			double ddp_dT_drho = _net->media->DpDTDrho(mid_rho, mid_T);
			double ddp_dT_dp = ddp_dT_drho / dp_drho;
			double dp_dT = _net->media->DpDT(mid_rho, mid_T);

			double B21 = cross_area - mid_m * mid_m / (cross_area * mid_rho * mid_rho) * drho_dp;
			double B22 = 2 * mid_m / cross_area / mid_rho;
			double F2 = -friction_coe / 2 * mid_m * fabs(mid_m) / (diameter * cross_area * mid_rho)
				- cross_area * mid_rho * GRAVITY * slope + mid_m * mid_m / cross_area / mid_rho / mid_rho * drho_dT * dT_dx;
			double G21 = mid_m * mid_m / cross_area * (2 / pow(mid_rho, 3)*drho_dp*drho_dp - 1 / mid_rho / mid_rho * ddrho_dp_dp)*dp_dx
				- 2 * mid_m / cross_area / mid_rho / mid_rho * drho_dp*dm_dx;
			double G22 = -2 * mid_m / cross_area / mid_rho / mid_rho * drho_dp*dp_dx
				+ 2 / cross_area / mid_rho * dm_dx;
			double dF_dU_21 = friction_coe / 2 * mid_m*fabs(mid_m) / diameter / cross_area / mid_rho / mid_rho * drho_dp - cross_area * GRAVITY * slope * drho_dp
				+ mid_m * mid_m / cross_area * dT_dx * (-2 / pow(mid_rho, 3) *drho_dp * drho_dT - 1 / mid_rho / mid_rho * (ddp_dT_dp / dp_drho + dp_dT * ddrho_dp_dp));
			double dF_dU_22 = -friction_coe * fabs(mid_m) / diameter / cross_area / mid_rho + 2 * mid_m / cross_area / mid_rho / mid_rho * drho_dT*dT_dx;

			double CE21 = 0 - 1 / dx * B21 + 0.5 * (G21 - dF_dU_21);
			double CE22 = 0.5 / dt - 1 / dx * B22 + 0.5*(G22 - dF_dU_22);
			double DW21 = 0 + 1 / dx * B21 + 0.5 * (G21 - dF_dU_21);
			double DW22 = 0.5 / dt + 1 / dx * B22 + 0.5*(G22 - dF_dU_22);
			double H2 = F2 + (G21 - dF_dU_21 + 0) * mid_p + (G22 - dF_dU_22 + 1 / dt) * mid_m;

			pipe_equation_coe_matrix(i, cur_node_index * 2) = CE21;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 1) = CE22;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 2) = DW21;
			pipe_equation_coe_matrix(i, cur_node_index * 2 + 3) = DW22;
			pipe_equation_cons_matrix(i, 0) = 0;
			pipe_equation_cons_matrix(i, 1) = 0;
			pipe_equation_cons_matrix(i, 2) = H2;
		}
	}
	//pipe_equation_coe_matrix.print("coe_matrix");
	//pipe_equation_cons_matrix.print("cons_matrix");

	//arma::mat temp(pipe_equation_coe_matrix);
	//for (int i = 0; i < matrix_order; i++) {
	//	for (int j = 0; j < matrix_order; j++) {
	//		if (temp(i, j) != 0) {
	//			printf("(%d, %d): %.5e\n", i, j, temp(i, j));
	//		}
	//	}
	//}

	return arma::spsolve(pipe_equation_coe_matrix, pipe_equation_cons_matrix);
}