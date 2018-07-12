#include "GasSource.h"

GasSource::GasSource(int comp_id, string source_name, int source_id, BoundaryConditionType bound_type, SourceDataQueue& p, SourceDataQueue& m, SourceDataQueue& T):
source_id(source_id), bound_type(bound_type), pressure(p), mass_flowrate(m), temperature(T)
{
	nonpipe_comp_id = comp_id;
	comp_name = source_name;
	comp_type = GasSource_Type;

	if (!pressure.empty() && pressure.top().time_occur > 0) {
		double value = pressure.top().value;
		pressure.push(SignalData("step", 0.0, value));
	}

	if (!mass_flowrate.empty() && mass_flowrate.top().time_occur > 0) {
		double value = mass_flowrate.top().value;
		mass_flowrate.push(SignalData("step", 0.0, value));
	}

	if (!temperature.empty() && temperature.top().time_occur > 0) {
		double value = temperature.top().value;
		temperature.push(SignalData("step", 0.0, value));
	}
}