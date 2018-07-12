#include "BWRS.h"

double pure_comp_coe_all[] = {
	200.7311407, 4.429759507462288e-2, 2171357.357, 36183311.07, 1288272176, 7.751494691, 5.278774066260021e-3, 297541.3871, 237.3927511, 6.890669426037629e-5, 5.360830749576498e-3,
	447.4268195, 6.740708589838079e-2, 18582413.44, 524514722.1, 16075898116, 31.18914617, 1.235845516730291e-2, 3312808.995, 2046.048796, 2.139682871434741e-4, 1.133369444348592e-2,
	701.1624142, 9.237418133819698e-2, 52566506.88, 1834141364, 52228091579, 74.12321148, 2.334412021890946e-2, 11981515.7, 6647.21561, 5.085427160949330e-4, 2.010275544550133e-2,
	996.1468918, 0.122281695, 99477895.29, 3858596945, 1.09686e+11, 146.0780207, 4.100898413510214e-2, 29185970.6, 15153.68351, 1.135920895422359e-3, 3.428616997037291e-2,
	994.3570248, 0.118948681, 112927137.3, 4579790905, 1.29236e+11, 145.4286983, 3.885343204893397e-2, 31777795.52, 16088.53656, 1.024393061068502e-3, 3.197148171211579e-2,
	1268.805286, 0.144686151, 185506562.5, 8200803414, 2.29278e+11, 237.611562, 5.763054745888307e-2, 61787856.45, 29757.73647, 1.766814760310705e-3, 4.588692452681199e-2,
	1277.753591, 0.147059985, 210599460.1, 9544688471, 2.54154e+11, 254.5718103, 5.966212331051422e-2, 69706592.38, 33704.36856, 1.785321042578454e-3, 4.612349655781080e-2,
	1562.310945, 0.176153268, 348506731.4, 17212990566, 4.44898e+11, 406.7452831, 8.591685033841574e-2, 132861063.9, 61821.82574, 2.847810283236869e-3, 6.275603892981481e-2,
	1836.409363, 0.206436062, 536118798.5, 28403010565, 7.18658e+11, 611.9521166, 0.118379889, 231178302.6, 104401.113, 4.242975754562693e-3, 8.159663009382195e-2,
	2080.139263, 0.238844329, 787834345.3, 44221029901, 1.09143e+12, 888.0619373, 0.158973249, 378935695.5, 167890.3797, 6.000446817656721e-3, 0.102441185,
	2273.084883, 0.270630404, 1108008720, 65417224498, 1.58298e+12, 1226.956639, 0.204674305, 583110492.7, 254171.4597, 7.904594978830571e-3, 0.122682756,
	2433.166781, 0.304636959, 1507644171, 92912399833, 2.20887e+12, 1658.419656, 0.259936245, 864473130.5, 371818.0598, 1.014998073655651e-2, 0.144493851,
	347.6328602, 5.645795589731836e-2, 12363925.75, 323356571.1, 9226757180, 20.26627382, 8.668908236424026e-3, 1848244.493, 1229.808853, 1.258654014430108e-4, 7.957354795836535e-3,
	629.6814305, 8.344024328668204e-2, 44843615.91, 1540745927, 44593064519, 59.37051411, 1.903378341497509e-2, 9306806.343, 5183.108661, 3.786217119562229e-4, 1.652200598879160e-2,
	118.3278445, 4.033973491124779e-2, 625766.9523, 7024090.547, 138123507.9, 4.345708704, 4.390476470497338e-3, 74572.14753, 96.13463608, 5.044216091363597e-5, 4.347002012372144e-3,
	259.3003579, 4.398705686483169e-2, 15716657.03, 457294750.1, 8853327136, 14.35154721, 5.319338926663501e-3, 1615134.22, 1159.124433, 5.082706286618077e-5, 4.313012393472809e-3,
	350.2941875, 4.330345289944436e-2, 22028690.86, 761323521.9, 28033610952, 15.78047787, 5.102151303404438e-3, 2511657.01, 1275.863507, 5.646859212462923e-5, 4.662007154237683e-3,
	291.5003376, 2.707086568756302e-2, 119055380.7, 7544594015, 2.31687e+11, 12.5416312, 2.034590209036903e-3, 6772430.275, 2540.144914, 9.700063586746393e-6, 1.416838415979742e-3
};

double interaction_coe[18][18] = {
	0.0000,0.0100,0.0230,0.0275,0.0310,0.0360,0.0410,0.0500,0.0600,0.0700,0.0810,0.0920,0.1010,0.0210,0.0250,0.0500,0.0500,0.0000,
	0.0100,0.0000,0.0031,0.0040,0.0045,0.0050,0.0060,0.0070,0.0085,0.0100,0.0120,0.0130,0.0150,0.0030,0.0700,0.0480,0.0450,0.0000,
	0.0230,0.0031,0.0000,0.0030,0.0035,0.0040,0.0045,0.0050,0.0065,0.0800,0.0100,0.0110,0.0130,0.0000,0.1000,0.0450,0.0400,0.0000,
	0.0275,0.0040,0.0030,0.0000,0.0000,0.0080,0.0010,0.0015,0.0018,0.0020,0.0025,0.0030,0.0030,0.0030,0.1100,0.0500,0.0360,0.0000,
	0.0310,0.0045,0.0035,0.0000,0.0000,0.0080,0.0010,0.0015,0.0018,0.0020,0.0025,0.0030,0.0030,0.0035,0.1200,0.0500,0.0340,0.0000,
	0.0360,0.0050,0.0040,0.0080,0.0080,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0040,0.1340,0.0500,0.0280,0.0000,
	0.0410,0.0060,0.0045,0.0010,0.0010,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0045,0.1480,0.0500,0.0200,0.0000,
	0.0500,0.0070,0.0050,0.0015,0.0015,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0050,0.1720,0.0500,0.0000,0.0000,
	0.0600,0.0085,0.0065,0.0018,0.0018,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0065,0.2000,0.0500,0.0000,0.0000,
	0.0700,0.0100,0.0800,0.0020,0.0020,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0080,0.2280,0.0500,0.0000,0.0000,
	0.0810,0.0120,0.0100,0.0025,0.0025,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0100,0.2640,0.0500,0.0000,0.0000,
	0.0920,0.0130,0.0110,0.0030,0.0030,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0110,0.2040,0.0500,0.0000,0.0000,
	0.0100,0.0000,0.0031,0.0040,0.0045,0.0050,0.0060,0.0070,0.0085,0.0100,0.0120,0.0130,0.0150,0.0030,0.0700,0.0480,0.0450,0.0000,
	0.0210,0.0030,0.0000,0.0030,0.0035,0.0040,0.0045,0.0050,0.0065,0.0080,0.0100,0.0110,0.0130,0.0000,0.1000,0.0450,0.0400,0.0000,
	0.0250,0.0700,0.1000,0.1100,0.1200,0.1340,0.1480,0.1720,0.2000,0.2280,0.2640,0.2040,0.3220,0.1000,0.0000,0.0000,0.0000,0.0000,
	0.0500,0.0480,0.0450,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0450,0.0000,0.0000,0.0350,0.0000,
	0.0500,0.0450,0.0400,0.0360,0.0340,0.0280,0.0200,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0400,0.0000,0.0350,0.0000,0.0000,
	0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000
};

double molecular_mass[] = { 16.042, 30.068, 44.094, 58.12, 58.12, 72.146, 72.146, 86.172, 100.198, 114.224, 128.25, 142.276, 28.05, 42.08, 28.016, 44.01, 34.076, 18.015 };

BWRS::BWRS(unordered_map<GasComponentType, double> comp){
	component = comp;
	equation_type = StateEquationType::BWRS;
	BWRS_coefficients.resize(11, 0);
	average_molecular_mass = 0;

	for (auto i = component.begin(); i != component.end(); i++) {
		BWRS_coefficients[1] += i->second * pure_comp_coe_all[i->first * 11 + 1];
		BWRS_coefficients[5] += i->second * pow(pure_comp_coe_all[i->first * 11 + 5], 1.0 / 3);
		BWRS_coefficients[6] += i->second * pow(pure_comp_coe_all[i->first * 11 + 6], 1.0 / 3);
		BWRS_coefficients[7] += i->second * pow(pure_comp_coe_all[i->first * 11 + 7], 1.0 / 3);
		BWRS_coefficients[8] += i->second * pow(pure_comp_coe_all[i->first * 11 + 8], 1.0 / 3);
		BWRS_coefficients[9] += i->second * pow(pure_comp_coe_all[i->first * 11 + 9], 1.0 / 3);
		BWRS_coefficients[10] += i->second * pow(pure_comp_coe_all[i->first * 11 + 10], 1.0 / 2);
		average_molecular_mass += molecular_mass[i->first] * i->second;

		for (auto j = component.begin(); j != component.end(); j++) {
			BWRS_coefficients[0] += i->second*j->second*sqrt(pure_comp_coe_all[i->first * 11])*sqrt(pure_comp_coe_all[j->first * 11])*(1 - interaction_coe[i->first][j->first]);
			BWRS_coefficients[2] += i->second*j->second*sqrt(pure_comp_coe_all[i->first * 11 + 2])*sqrt(pure_comp_coe_all[j->first * 11 + 2]) * pow(1 - interaction_coe[i->first][j->first], 3);
			BWRS_coefficients[3] += i->second*j->second*sqrt(pure_comp_coe_all[i->first * 11 + 3])*sqrt(pure_comp_coe_all[j->first * 11 + 3]) * pow(1 - interaction_coe[i->first][j->first], 4);
			BWRS_coefficients[4] += i->second*j->second*sqrt(pure_comp_coe_all[i->first * 11 + 4])*sqrt(pure_comp_coe_all[j->first * 11 + 4]) * pow(1 - interaction_coe[i->first][j->first], 5);
		}
	}
	BWRS_coefficients[5] = pow(BWRS_coefficients[5], 3);
	BWRS_coefficients[6] = pow(BWRS_coefficients[6], 3);
	BWRS_coefficients[7] = pow(BWRS_coefficients[7], 3);
	BWRS_coefficients[8] = pow(BWRS_coefficients[8], 3);
	BWRS_coefficients[9] = pow(BWRS_coefficients[9], 3);
	BWRS_coefficients[10] = pow(BWRS_coefficients[10], 2);

}

double BWRS::CalculateDensity(double p, double T) {
	p /= 1000;//ѹ������λת��KPa, ���ڼ���
	double rho1 = 10, rho2 = p / R_CONSTANT / T;

	double T_2 = T * T;
	double T_3 = T * T * T;
	double T_4 = T_2 * T_2;
	
	while (fabs(rho1 - rho2) > 1e-6) {
		double rho1_2 = rho1 * rho1;
		double rho1_3 = rho1 * rho1 * rho1;
		double rho1_6 = rho1_3 * rho1_3;
		double rho2_2 = rho2 * rho2;
		double rho2_3 = rho2 * rho2 * rho2;
		double rho2_6 = rho2_3 * rho2_3;

		double F1 = rho1 * R_CONSTANT*T + (BWRS_coefficients[1] * R_CONSTANT*T - BWRS_coefficients[0] - BWRS_coefficients[2] / T_2 + BWRS_coefficients[3] / T_3 - BWRS_coefficients[4] / T_4)*rho1_2
			+ (BWRS_coefficients[6] * R_CONSTANT * T - BWRS_coefficients[5] - BWRS_coefficients[8] / T) * rho1_3
			+ BWRS_coefficients[9] * (BWRS_coefficients[5] + BWRS_coefficients[8] / T) * rho1_6
			+ BWRS_coefficients[7] * rho1_3 / T_2 * (1 + BWRS_coefficients[10] * rho1_2) * exp(-BWRS_coefficients[10] * rho1_2) - p;
		double F2 = rho2 * R_CONSTANT*T + (BWRS_coefficients[1] * R_CONSTANT*T - BWRS_coefficients[0] - BWRS_coefficients[2] / T_2 + BWRS_coefficients[3] / T_3 - BWRS_coefficients[4] / T_4)*rho2_2
			+ (BWRS_coefficients[6] * R_CONSTANT * T - BWRS_coefficients[5] - BWRS_coefficients[8] / T) * rho2_3
			+ BWRS_coefficients[9] * (BWRS_coefficients[5] + BWRS_coefficients[8] / T) * rho2_6
			+ BWRS_coefficients[7] * rho2_3 / T_2 * (1 + BWRS_coefficients[10] * rho2_2) * exp(-BWRS_coefficients[10] * rho2_2) - p;

		double temp = rho2;
		rho2 = (rho1 * F2 - rho2 * F1) / (F2 - F1);
		rho1 = temp;
	}

	return rho2 * average_molecular_mass;
}

double BWRS::DpDrho(double rho, double T) {
	rho = rho / average_molecular_mass;  ////ת����λKmol / m3
	double T_2 = T * T;
	double T_3 = T_2 * T;
	double T_4 = T_2 * T_2;
	double rho_2 = rho * rho;
	double rho_3 = rho * rho_2;
	double rho_4 = rho_2 * rho_2;
	double result = R_CONSTANT * T + 2.0*(BWRS_coefficients[1] * R_CONSTANT * T - BWRS_coefficients[0] - BWRS_coefficients[2] / T_2 +
		BWRS_coefficients[3] / T_3 - BWRS_coefficients[4] / T_4)*rho +
		3.0*(BWRS_coefficients[6] * R_CONSTANT * T - BWRS_coefficients[5] - BWRS_coefficients[8] / T)*rho_2 +
		6.0*BWRS_coefficients[9] * (BWRS_coefficients[5] + BWRS_coefficients[8] / T)*rho_2*rho_3 +
		3.0*BWRS_coefficients[7] * rho_2 / T_2 * (1.0 + BWRS_coefficients[10] * rho_2 -
			2.0 / 3.0*BWRS_coefficients[10] * BWRS_coefficients[10] * rho_4)*
		exp(-BWRS_coefficients[10] * rho_2);   /////��BWRS���̼���DPDRHO�ĵ�λΪKPa / (Kmol / m3)
	return result * 1000.0 / average_molecular_mass;  ////ת����λPa / (Kg / m3)
}

double BWRS::DpDrhoDrho(double rho, double T) {
	rho = rho / average_molecular_mass;//ת����λKmol / m3

	double T_2 = T * T;
	double T_3 = T_2 * T;
	double T_4 = T_2 * T_2;
	double rho_2 = rho * rho;
	double rho_3 = rho * rho_2;
	double rho_6 = rho_3 * rho_3;

	double result = T * (2.0*BWRS_coefficients[1] * R_CONSTANT + 6.0*rho*BWRS_coefficients[6] * R_CONSTANT)
		- 2.0*BWRS_coefficients[0] + (2.0*BWRS_coefficients[3]) / T_3 - (2.0*BWRS_coefficients[4]) / T_4
		- 6.0*rho*BWRS_coefficients[5] + (2.0*(3.0*rho*BWRS_coefficients[7] - BWRS_coefficients[2]
			* exp(rho_2*BWRS_coefficients[10]) + 3.0*rho_3*BWRS_coefficients[7] * BWRS_coefficients[10]
			- 9.0*rho_2*rho_3*BWRS_coefficients[7] * BWRS_coefficients[10] * BWRS_coefficients[10] + 2.0*rho_6*rho*BWRS_coefficients[7]
			* BWRS_coefficients[10] * BWRS_coefficients[10] * BWRS_coefficients[10])) / (T_2*exp(rho_2*BWRS_coefficients[10]))
		+ 30.0*rho_2*rho_2*BWRS_coefficients[5] * BWRS_coefficients[9] + (6.0*rho*BWRS_coefficients[8]
			* (5.0*rho_3*BWRS_coefficients[9] - 1.0)) / T;
	/////��BWRS���̼���DPDRHODRHO�ĵ�λΪKPa / (Kmol / m3) / (Kmol / m3)
	return result * 1000.0 / average_molecular_mass / average_molecular_mass; ////ת����λPa / (Kg / m3) / (Kg / m3)
}

double BWRS::DpDT(double rho, double T) {

	rho = rho / average_molecular_mass;  ////ת����λKmol / m3
	double T_2 = T * T;
	double T_3 = T_2 * T;
	double T_4 = T_2 * T_2;
	double rho_2 = rho * rho;
	double rho_3 = rho * rho_2;
	double rho_6 = rho_3 * rho_3;
	double result = rho * R_CONSTANT + (BWRS_coefficients[1] * R_CONSTANT + 2.0*BWRS_coefficients[2] / T_3
		- 3.0*BWRS_coefficients[3] / T_4 + 4.0*BWRS_coefficients[4] / T_4 / T)*rho_2 +
		(BWRS_coefficients[6] * R_CONSTANT + BWRS_coefficients[8] / T_2)*rho_3 -
		BWRS_coefficients[9] * (BWRS_coefficients[8] / T_2)*rho_6 -
		2.0*BWRS_coefficients[7] * rho_3 / T_3 * (1.0 + BWRS_coefficients[10] * rho_2)*
		exp(-BWRS_coefficients[10] * rho_2);
	/////��BWRS���̼���DPDT�ĵ�λΪKPa / K
	return result * 1000.0;   ////ת����λPa / K
}

double BWRS::DpDTDT(double rho, double T) {
	rho = rho / average_molecular_mass;  ////ת����λKmol / m3
	double T_2 = T * T;
	double T_3 = T_2 * T;
	double T_4 = T_2 * T_2;
	double rho_2 = rho * rho;
	double rho_3 = rho * rho_2;
	double rho_6 = rho_3 * rho_3;
	double result = (-6.0*BWRS_coefficients[2] / T_4 +
		12.0*BWRS_coefficients[3] / T_3 / T_2 - 20.0*BWRS_coefficients[4] / T_3 / T_3)*rho_2 +
		(-2.0*BWRS_coefficients[8] / T_3)*rho_3 +
		2.0*BWRS_coefficients[9] * (BWRS_coefficients[8] / T_3)*rho_6 +
		6.0*BWRS_coefficients[7] * rho_3 / T_4 * (1.0 + BWRS_coefficients[10] * rho_2)*
		exp(-BWRS_coefficients[10] * rho_2);
	/////��BWRS���̼���DPDT�ĵ�λΪKPa / K / K
	return result * 1000.0;   ////ת����λPa / K / K
}

double BWRS::DpDTDrho(double rho, double T) {
	rho = rho / average_molecular_mass;  ////ת����λKmol / m3
	double T_2 = T * T;
	double T_3 = T_2 * T;
	double rho_2 = rho * rho;
	double rho_3 = rho * rho_2;
	double rho_5 = rho_2 * rho_3;
	double result = R_CONSTANT + 2.0*rho*BWRS_coefficients[1] * R_CONSTANT - (6.0*rho*BWRS_coefficients[3]) / T_2 / T_2
		+ (8.0*rho*BWRS_coefficients[4]) / T_2 / T_3 + 3.0*rho_2*BWRS_coefficients[6] * R_CONSTANT
		- (3.0*rho_2*BWRS_coefficients[8] * (2.0*rho_3*BWRS_coefficients[9] - 1.0)) / T_2
		+ (2.0*rho*(2.0*BWRS_coefficients[2] * exp(rho_2*BWRS_coefficients[10]) - 3.0*rho*BWRS_coefficients[7]
			- 3.0*rho_3*BWRS_coefficients[7] * BWRS_coefficients[10] + 2.0*rho_5*BWRS_coefficients[7] * BWRS_coefficients[10] * BWRS_coefficients[10]))
		/ (T_3*exp(rho_2*BWRS_coefficients[10]));
	/////��BWRS_coefficients���̼���DPDTDrho�ĵ�λΪKPa / K / (Kmol / m3)
	return result * 1000.0 / average_molecular_mass;    ////ת����λPa / K / (Kg / m3)
}

double BWRS::DrhoDT(double rho, double T) {
	double dp_dT = this->DpDT(rho, T);
	double dp_drho = this->DpDrho(rho, T);
	return -dp_dT / dp_drho;
}