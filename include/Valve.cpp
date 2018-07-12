#include "Valve.h"

Valve::Valve() {

}

Valve::Valve(int comp_id, string valve_name, int valve_id, ValveDataQueue valve_status) :
	valve_id(valve_id), status(valve_status) {
	nonpipe_comp_id = comp_id;
	comp_name = valve_name;
	is_reverse = false;

	if (status.top().time_occur > 0) {
		string status_mode = status.top().signal_mode;
		
		if (status_mode == "on") {
			double pressure_ratio = status.top().value;
			status.push(SignalData(status_mode, 0, pressure_ratio));
		}
		else {
			status.push(SignalData(status_mode, 0, 0));
		}
	}
}