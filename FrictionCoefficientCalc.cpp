#include "FrictionCoefficientCalc.h"

double FrictionCoefficientCalc(Pipe* cur_pipe, double rho, double v, double dy_vis) {
	if (cur_pipe->wall_roughness > 1e-6) {
		double inter_diameter = cur_pipe->inter_diameter;
		double Re = rho * fabs(v) * inter_diameter / dy_vis;

		if (Re <= 3000.0)
			return 0.1 - 2.442857142857143e-5*Re;
		else {
			double res1 = 100;
			double res2 = 7;

			while (fabs(res1 - res2) > 1e-6) {
				double F1 = res1 + 2 * log10(cur_pipe->wall_roughness / 3.7 / inter_diameter + 2.51 / Re * res1);
				double F2 = res2 + 2 * log10(cur_pipe->wall_roughness / 3.7 / inter_diameter + 2.51 / Re * res2);

				double temp = res2;
				res2 = (F2*res1 - F1 * res2) / (F2 - F1);
				res1 = temp;
			}
			return 1.0 / res2 / res2;
		}
	}
	else
		return 0.0;
}