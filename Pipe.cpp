#include "Pipe.h"
#include <algorithm>

#define PI 3.1415926535898

Pipe::Pipe(int id, string name, int inlet_id, int outlet_id, double len, double outerD, 
	double thickness, double roughness, vector<ElevationData> ele, int node_num, 
	vector<InitialPipeData> init_p, vector<InitialPipeData> init_m, vector<InitialPipeData> init_T) :
	pipe_id(id),
	pipe_name(name),
	inlet_comp_id(inlet_id),
	outlet_comp_id(outlet_id),
	length(len),
	outer_diameter(outerD),
	wall_thickness(thickness),
	wall_roughness(roughness),
	elevation_set(ele),
	node_num(node_num)
{
	inter_diameter = outer_diameter - wall_thickness * 2;
	cross_area = PI * inter_diameter * inter_diameter / 4;
	dx = length / (node_num - 1);

	sort(elevation_set.begin(), elevation_set.end());
	auto new_tail_ele = unique(elevation_set.begin(), elevation_set.end());
	elevation_set.erase(new_tail_ele, elevation_set.end());

	sort(init_p.begin(), init_p.end());
	auto new_tail_p = unique(init_p.begin(), init_p.end());
	init_p.erase(new_tail_p, init_p.end());

	sort(init_m.begin(), init_m.end());
	auto new_tail_m = unique(init_m.begin(), init_m.end());
	init_m.erase(new_tail_m, init_m.end());

	sort(init_T.begin(), init_T.end());
	auto new_tail_T = unique(init_T.begin(), init_T.end());
	init_T.erase(new_tail_T, init_T.end());

	vector<double> p(node_num, 0);
	vector<double> m(node_num, 0);
	vector<double> T(node_num, 0);
	for (int i = 0; i < node_num; i++) {
		double horizontal_distance = i * dx;

		auto iter_high_p = upper_bound(init_p.begin(), init_p.end(), horizontal_distance);
		auto iter_high_m = upper_bound(init_m.begin(), init_m.end(), horizontal_distance);
		auto iter_high_T = upper_bound(init_T.begin(), init_T.end(), horizontal_distance);

		if (iter_high_p == init_p.end()) p[i] = (iter_high_p - 1)->value;
		else if (iter_high_p - 1 == init_p.begin()) p[i] = iter_high_p->value;
		else {
			double segment_len = iter_high_p->horizontal_distance - (iter_high_p - 1)->horizontal_distance;
			double pressure_difference = iter_high_p->value - (iter_high_p - 1)->value;
			double change_rate = pressure_difference / segment_len;
			p[i] = (iter_high_p - 1)->value + change_rate * (horizontal_distance - (iter_high_p - 1)->horizontal_distance);
		}

		if (iter_high_m == init_m.end()) m[i] = (iter_high_m - 1)->value;
		else if (iter_high_m - 1 == init_m.begin()) m[i] = iter_high_m->value;
		else {
			double segment_len = iter_high_m->horizontal_distance - (iter_high_m - 1)->horizontal_distance;
			double massflowrate_difference = iter_high_m->value - (iter_high_m - 1)->value;
			double change_rate = massflowrate_difference / segment_len;
			m[i] = (iter_high_m - 1)->value + change_rate * (horizontal_distance - (iter_high_m - 1)->horizontal_distance);
		}

		if (iter_high_T == init_T.end()) T[i] = (iter_high_T - 1)->value;
		else if (iter_high_T - 1 == init_T.begin()) T[i] = iter_high_T->value;
		else {
			double segment_len = iter_high_T->horizontal_distance - (iter_high_T - 1)->horizontal_distance;
			double temperature_difference = iter_high_T->value - (iter_high_T - 1)->value;
			double change_rate = temperature_difference / segment_len;
			T[i] = (iter_high_T - 1)->value + change_rate * (horizontal_distance - (iter_high_T - 1)->horizontal_distance);
		}
	}

	pressure = p;
	massflowrate = m;
	temperature.push_back(T);
	temperature.push_back(T);
}

double Pipe::GetSlope(int node_index) {
	double horizontal_distance = node_index * dx;

	auto iter_high = upper_bound(elevation_set.begin(), elevation_set.end(), horizontal_distance);
	if (iter_high == elevation_set.end() || iter_high - 1 == elevation_set.begin()) {
		return 0;
	}
	double elevation_difference = (iter_high->elevation - (iter_high - 1)->elevation);
	double segment_len = iter_high->horizontal_distance - (iter_high - 1)->horizontal_distance;

	double slope = elevation_difference / sqrt(segment_len * segment_len + elevation_difference * elevation_difference);

	return slope;
}