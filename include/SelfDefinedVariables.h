#pragma once

#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

enum NonPipeComponentType {
	Valve_Type,
	Compressor_Type,
	GasSource_Type
};

enum BoundaryConditionType {
	None,
	OnlyPressure,
	OnlyMassFlowRate,
	BothPressureAndMassFlowRate,
};

enum GasComponentType
{
	CH4,
	C2H6,
	C3H8,
	i_C4H10,
	n_C4H10,
	i_C5H12,
	n_C5H12,
	C6H14,
	C7H16,
	C8H18,
	C9H20,
	C10H22,
	C2H4,
	C3H6,
	N2,
	CO2,
	H2S,
	H2O,
};

struct SignalData {
	string signal_mode;
	double time_occur;
	double value;

	SignalData() {}
	SignalData(string s, double t, double v) :signal_mode(s), time_occur(t), value(v) {}

	friend bool operator > (const SignalData& a, const SignalData& b) {
		return a.time_occur > b.time_occur;
	}
	friend bool operator < (const SignalData& a, const SignalData& b) {
		return a.time_occur < b.time_occur;
	}
	friend bool operator == (const SignalData& a, const SignalData& b) {
		return a.time_occur == b.time_occur;
	}
};

struct ElevationData {
	double horizontal_distance;
	double elevation;

	ElevationData(){}
	ElevationData(double h, double e):horizontal_distance(h),elevation(e){}

	friend bool operator < (const ElevationData& a, const ElevationData& b) {
		return a.horizontal_distance < b.horizontal_distance;
	}
	friend bool operator > (const ElevationData& a, const ElevationData& b) {
		return a.horizontal_distance > b.horizontal_distance;
	}
	friend bool operator == (const ElevationData& a, const ElevationData& b) {
		return a.horizontal_distance == b.horizontal_distance;
	}
	friend bool operator < (const double& b, const ElevationData& a) {
		return b < a.horizontal_distance;
	}
	friend bool operator > (const double& b, const ElevationData& a) {
		return b > a.horizontal_distance;
	}
};

struct InitialPipeData {
	double horizontal_distance;
	double value;

	InitialPipeData(){}
	InitialPipeData(double h,double v):horizontal_distance(h),value(v){}

	friend bool operator < (const InitialPipeData& a, const InitialPipeData& b) {
		return a.horizontal_distance < b.horizontal_distance;
	}
	friend bool operator > (const InitialPipeData& a, const InitialPipeData& b) {
		return a.horizontal_distance > b.horizontal_distance;
	}
	friend bool operator == (const InitialPipeData& a, const InitialPipeData& b) {
		return a.horizontal_distance == b.horizontal_distance;
	}
	friend bool operator < (const double& b, const InitialPipeData& a) {
		return b < a.horizontal_distance;
	}
	friend bool operator > (const double& b, const InitialPipeData& a) {
		return b > a.horizontal_distance;
	}
};