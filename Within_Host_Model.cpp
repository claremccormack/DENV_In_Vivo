#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include<random>
#include<algorithm>
#include<time.h>
#include<fstream>
#include<string>
#include<sstream>
#include<direct.h>
#include<windows.h>
#include<limits.h>
#include <omp.h>
#include"randlib_par.h"
#include"Structures.h"
#include"Functions.h"
#include"Within_Host_Model.h"

using namespace std;

void Within_Host_Model(host &Host, int thread, int start, int steps, int runs, model_type model) {

		double
			h_x,h_y, h_v,h_m,h_g,h_n,
			p_x,p_y, p_v,p_m,p_g,p_n,
			t1;

		long double
			current_x, current_y, current_v, current_m, current_g, current_d, current_n,
			new_x, new_y, new_v, innoc_v, new_n, new_m, new_g,
			out_x,out_y, out_v, out_n, out_m, out_g,
			deaths_x, deaths_y, deaths_v, cleared_y, cleared_v;

		int m = thread;
	
		model_params params = Host.Params;

		for (int r = 0; r < runs; r++) {

			for (int t = start; t < steps + 1; t++)

			{
				t1 = int(t / (1 / DT)); /*day number*/

				/*current population in each compartment*/

				current_x = Host.X[m][t][r];
				current_y = Host.Y[m][t][r];
				current_v = Host.V[m][t][r];
				current_n = Host.N[m][t][r];
				current_m = Host.M[m][t][r];
				current_g = Host.G[m][t][r]; 
				current_d = Host.D[m][t][r];

				/*hazard of leaving each compartment*/

					h_x = params.gamma +Infection_Rate(Host,thread,t,r,model,params.beta,Host.PK.emax,Host.Params.ec50,Host.Params.hill_coeff);
					h_y = params.delta + params.alpha * current_m * (1 / (1 + params.c * current_m));
				
					h_v = params.kappa;
					h_n = params.kappa_n;
					h_m = params.kappa_m;
					h_g = params.kappa_g; 
					
				/*convert hazards to probabilities*/

					p_x = 1 - exp(-h_x * DT);
					p_y = 1 - exp(-h_y * DT);
					p_v = 1 - exp(-h_v * DT);
					p_n = 1 - exp(-h_n * DT);
					p_m = 1 - exp(-h_m * DT);
					p_g = 1 - exp(-h_g * DT); 
					
				/*determine movement between compartments*/

						new_x = params.A * DT;
						out_x = current_x * p_x;
						deaths_x = out_x * params.gamma / h_x;
						if (params.gamma == h_x) { new_y = 0; }
						else { new_y = out_x - deaths_x; }

						out_y = current_y * p_y;
						deaths_y = out_y * (params.delta / h_y);
						cleared_y = out_y - deaths_y;

						new_v = Virus_Proliferation(Host, thread, t, r, model, params.omega,Host.PK.emax, Host.Params.ec50, Host.Params.hill_coeff);
						out_v = current_v * p_v;
						deaths_v = out_v * (params.kappa / h_v);
						cleared_v = out_v - deaths_v;

						if (t == Host.time_infection + 24 * RUN_IN && t > 0) { innoc_v = Host.scaling_inoculum * Host.viral_inoculum; }

						else { innoc_v = 0; }

						new_n = Antigen_Proliferation(Host, thread, t, r, model, params.rho,Host.PK.emax, Host.Params.ec50, Host.Params.hill_coeff);
						out_n = current_n * p_n;

						new_m = IgM_Response(Host, thread, t, r, model, params.eta,params.eta2);
						out_m = current_m * p_m;
										
						new_g= IgG_Response(Host, thread, t, r, model, params.psi,params.psi2);
						out_g = current_g * p_g;

				/*calculate population size at next timestep*/

				Host.X[m][t + 1][r] = current_x + new_x - deaths_x - new_y; 
				Host.Y[m][t + 1][r] = current_y + new_y - deaths_y - cleared_y;
				Host.V[m][t + 1][r] = current_v + innoc_v + new_v - deaths_v - cleared_v;
				Host.N[m][t + 1][r] = current_n + new_n - out_n;
				Host.M[m][t + 1][r] = current_m + new_m - out_m;
				Host.G[m][t + 1][r] = current_g + new_g - out_g;
				Host.D[m][t + 1][r] = PK_Model(Host, m, t, r, model);
			}
	}
}

double Infection_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	double new_rate, conc, effect;
	int m = thread;

	model_params params = Host.Params;

	conc = Host.D[m][t][r];

	if (model.drug_moa == 0) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

	else { effect = 0; }

	new_rate = Host.V[m][t][r]*rate*(1 - effect);

	return new_rate;
}

long double Virus_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	long double new_v, conc, effect,prod_rate;
	int m = thread;
	int t1;

	model_params params = Host.Params;

	conc = Host.D[m][t][r];

	if (model.drug_moa == 1) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

	else { effect = 0; }

	t1 = t;

	if (Host.V[m][t1][r] < model.threshold) { prod_rate = 0; }

	else { prod_rate = rate; }

	new_v = prod_rate * (1 - effect) * Host.Y[m][t1][r] * DT;

	return new_v;
}

long double Antigen_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	long double new_n, effect, conc;
	int m = thread;
	int t1;

	model_params params = Host.Params;

	new_n = 0;

	t1 = t;

	conc = Host.D[m][t][r];

	if (model.drug_moa == 1) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

	else { effect = 0; }

	new_n = rate *(1-effect)* Host.Y[m][t1][r] * DT;
		
	return new_n;

}

long double IgM_Response(host Host, int thread, int t, int r, model_type model, double p1, double p2) {

	long double new_m;
	int m = thread;
	int t1;

	model_params params = Host.Params;

	new_m = 0;

	t1 = t;

	new_m += (p1*Host.Y[m][t1][r]*Host.M[m][t][r] * DT) / (p2+ Host.Y[m][t1][r]);
		
	return new_m;

}

long double IgG_Response(host Host, int thread, int t, int r, model_type model, double p1, double p2) {

	long double new_g;
	int m = thread;
	int t1;

	model_params params = Host.Params;

	if (t < Host.Params.lag_g * 24) { t1 = 0; }

	else { t1 = t - round(Host.Params.lag_g * 24); }

	new_g = (p1 * Host.Y[m][t1][r] * Host.G[m][t][r] * DT) / (p2 + Host.Y[m][t1][r]);
		
	return new_g;

}

double PK_Model(host &Host, int thread, int t, int r, model_type model) {

	int m = thread;

	double cl_effect,current_A1,c_ss,clearance_A2,eliminated_A1;

	switch (model.pk_model_number) {

	case 1:

		/*AUC model (NHP)*/

		current_A1 = Host.PK.A1[m][t][r];

		if (t == Host.first_dose + 24 * RUN_IN) { c_ss = Host.PK.auc / 24; }

		else { c_ss = 0;}

		if (t > Host.last_dose+24* RUN_IN + 24) { eliminated_A1 = Host.PK.k_e*current_A1; }

		else { eliminated_A1 = 0; }

		Host.PK.A1[m][t + 1][r] = current_A1 + c_ss - eliminated_A1;

		return Host.PK.A1[m][t + 1][r];

		break;


	case 2:

		/*Janssen PK Model (Mice)*/

		if (Host.dose_frequency == 1) {

			cl_effect = (Host.PK.cl_emax1) / (1 + pow((Host.PK.cl_ed50 / Host.dose_size), 1));

			clearance_A2 = Host.PK.cl_max * (1 - cl_effect) * exp(Host.PK.cl_var);


		}

		else {

			cl_effect = (Host.PK.cl_emax2) / (1 + pow((Host.PK.cl_ed50 / Host.dose_size), 1));

			clearance_A2 = Host.PK.cl_max * (1 - 0.225) * (1 - cl_effect) * exp(Host.PK.cl_var);
		}

		current_A1 = Host.PK.A1[m][t][r];

		Host.PK.auc = Host.dose_size / clearance_A2;

		Host.PK.k_e = 1 - exp(-clearance_A2 / Host.PK.V);

		if (t == Host.first_dose + 24 * RUN_IN) {

			if (Host.dose_frequency == 1) { c_ss = Host.PK.auc * 1e6 / 24; }

			else { c_ss = Host.PK.auc * 1e6 / 12; }

		}

		else { c_ss = 0; }

		if (t > Host.last_dose + 24 * RUN_IN + 24) { eliminated_A1 = Host.PK.k_e * current_A1; }

		else { eliminated_A1 = 0; }

		Host.PK.A1[m][t + 1][r] = current_A1 + c_ss - eliminated_A1;

		return Host.PK.A1[m][t + 1][r];

		break;

	default:

		return 0;
	}

}

double Drug_Absorption(host &Host, int t) {

	double result;

	result = Host.dosage[t];

	return result;
}
