#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<cmath>
#include<ctime>
#include<random>
#include<algorithm>
#include<string>
#include<fstream>
#include"randlib_par.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_math.h"
#include<gsl/gsl_statistics.h>
#include"Structures.h"
#include"Functions.h"
#include"Within_Host_Model.h"
#include"MCMC.h"

random_device rd1;
mt19937 gen1(rd1());

gsl_ran_discrete_t * F1;
gsl_rng * g2 = gsl_rng_alloc(gsl_rng_mt19937);

const double na = sqrt(-1);

void Read_Data(string filename, expt*& Expt, int& n_expt) {

	ifstream Data_File;

	Data_File.open(filename);

	string line, experiment, group, individual, id, dose, frequency, inoculum, time_infection, time_measurement, v,  ns1, d, IgM, IgG, lod_v,  lod_m, lod_g, lod_n, time_death, cl_var, auc,species;

	int i, j, k, prev_i, prev_j, prev_k;

	int t = 0;

	prev_i = 0;
	prev_j = 0;
	prev_k = 0;

	n_expt = 0;

	/*count how many experiments*/

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);

			if (i == 0) { n_expt = 1; }

			if (i > prev_i) {

				n_expt += 1;
			}

			prev_i = i;
		}

		Data_File.close();

		Expt = new expt[n_expt];


	/*count how many groups within each experiment*/

		prev_i = 0;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);

			if (j == 0) { Expt[i].no_groups = 1; }

			if (i == prev_i && j > prev_j) {

				Expt[i].no_groups += 1;
			}

			prev_i = i;
			prev_j = j;
		}

		Data_File.close();

		for (int i = 0; i < n_expt; i++) { Expt[i].Pop = new pop[Expt[i].no_groups];}
	
	/*count how many individuals within each group*/

		prev_i = 0;
		prev_j = 0;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\t');
			getline(Data_File, individual, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);
			k = stoi(individual);

			if (k == 0) { Expt[i].Pop[j].no_hosts = 1; }

			if (j == prev_j && k > prev_k) {

				Expt[i].Pop[j].no_hosts += 1;
			}

			prev_i = i;
			prev_j = j;
			prev_k = k;
		}

		Data_File.close();

	/*allocate memory and sort data for each individual*/

		for (int i = 0; i < n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				Expt[i].Pop[j].Hosts = new host[Expt[i].Pop[j].no_hosts];
				Expt[i].Pop[j].Hosts_Avg = new host_avg[Expt[i].Pop[j].no_hosts];
			}
		}


		Expt[0].Pop[0].Hosts[0].Data.no_measurement_times = 1;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\t');
			getline(Data_File, individual, '\t');
			getline(Data_File, id, '\t');
			getline(Data_File, dose, '\t');
			getline(Data_File, frequency, '\t');
			getline(Data_File, inoculum, '\t');
			getline(Data_File, time_measurement, '\t');
			getline(Data_File, time_infection, '\t');
			getline(Data_File, v, '\t');
			getline(Data_File, ns1, '\t');
			getline(Data_File, IgM, '\t');
			getline(Data_File, IgG, '\t');
			getline(Data_File, d, '\t');
			getline(Data_File, lod_v, '\t');
			getline(Data_File, lod_n, '\t'); 
			getline(Data_File, lod_m, '\t');
			getline(Data_File, lod_g, '\t');
			getline(Data_File, cl_var, '\t');
			getline(Data_File, auc, '\t');
			getline(Data_File, species, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);
			k = stoi(individual);

			if (i != prev_i || j != prev_j || k != prev_k) {

				t = 0;

				Expt[i].Pop[j].Hosts[k].Data.no_measurement_times = 1;
			}

			if (t > 0 && i == prev_i && j == prev_j && k == prev_k) {

				Expt[i].Pop[j].Hosts[k].Data.no_measurement_times += 1;
			}

			Expt[i].Pop[j].Hosts[k].id = id;
			Expt[i].Pop[j].Hosts[k].dose_size = stod(dose);
			Expt[i].Pop[j].Hosts[k].dose_frequency = stod(frequency);
			Expt[i].Pop[j].Hosts[k].time_infection = stoi(time_infection);
			Expt[i].Pop[j].Hosts[k].viral_inoculum = stod(inoculum);
			Expt[i].Pop[j].Hosts[k].Data.lod_v = pow(10, stod(lod_v));/*limit of detection (extracellular virus)*/
			
			/*limit of detection (IgM)*/
			if (lod_m.compare("NA") == 0) {Expt[i].Pop[j].Hosts[k].Data.lod_m = 0;}
			else { Expt[i].Pop[j].Hosts[k].Data.lod_m = stod(lod_m);}

			/*limit of detection (NS1)*/
			if (lod_n.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.lod_n = 0; }
			else { Expt[i].Pop[j].Hosts[k].Data.lod_n = stod(lod_n); }

			/*limit of detection (IgG)*/
			if (lod_g.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.lod_g = 0; }
			else { Expt[i].Pop[j].Hosts[k].Data.lod_g = stod(lod_g); }

			Expt[i].Pop[j].Hosts[k].Data.Time[t] = stoi(time_measurement);

			if (v.compare("NA") == 0) {Expt[i].Pop[j].Hosts[k].Data.V[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.V[t] = stod(v);}

			if (d.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.D[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.D[t] = stod(d); }

			if (IgM.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.M[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.M[t] = stod(IgM); }

			if (IgG.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.G[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.G[t] = stod(IgG); }

			if (ns1.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.N[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.N[t] = stod(ns1); }

			if (cl_var.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].PK.cl_var = 0.0743; }
			else { Expt[i].Pop[j].Hosts[k].PK.cl_var = stod(cl_var); }

			if (species.compare("NHP") == 0) { Expt[i].Pop[j].Hosts[k].PK.cl_var = 0; }

			if (auc.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].PK.auc = 0; }
			else { Expt[i].Pop[j].Hosts[k].PK.auc = stod(auc);}

			prev_i = i;
			prev_j = j;
			prev_k = k;

			t = t + 1;
		}

	Data_File.close();
}

void Read_Param_Values(string filename, model_type& model, fit_params& Fit, expt*Expt) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				ifstream Param_File;

				Param_File.open(filename);

				string param_name, sample, param_value_string;

				while (!Param_File.eof())
				{
					getline(Param_File, param_name, '\t');
					getline(Param_File, param_value_string, '\n');

						if (param_name == "n_expt")					model.n_expt = stoi(param_value_string);
						
						if (param_name == "A")						Expt[i].Pop[j].Hosts[k].Params.A = stod(param_value_string);
						if (param_name == "alpha")					Expt[i].Pop[j].Hosts[k].Params.alpha= stod(param_value_string);
						if (param_name == "beta")					Expt[i].Pop[j].Hosts[k].Params.beta = stod(param_value_string);
						if (param_name == "delta")					Expt[i].Pop[j].Hosts[k].Params.delta = stod(param_value_string);
						if (param_name == "eta")					Expt[i].Pop[j].Hosts[k].Params.eta = stod(param_value_string);
						if (param_name == "eta2")					Expt[i].Pop[j].Hosts[k].Params.eta2 = stod(param_value_string);
						if (param_name == "gamma")					Expt[i].Pop[j].Hosts[k].Params.gamma = stod(param_value_string);
						if (param_name == "kappa")					Expt[i].Pop[j].Hosts[k].Params.kappa = stod(param_value_string);
						if (param_name == "kappa_m")				Expt[i].Pop[j].Hosts[k].Params.kappa_m = stod(param_value_string);
						if (param_name == "kappa_g")				Expt[i].Pop[j].Hosts[k].Params.kappa_g = stod(param_value_string);
						if (param_name == "kappa_n")				Expt[i].Pop[j].Hosts[k].Params.kappa_n = stod(param_value_string);
						if (param_name == "omega")					Expt[i].Pop[j].Hosts[k].Params.omega = stod(param_value_string);
						if (param_name == "rho")					Expt[i].Pop[j].Hosts[k].Params.rho = stod(param_value_string);
						if (param_name == "scaling_inoculum")		Expt[i].Pop[j].Hosts[k].scaling_inoculum = stod(param_value_string);
						if (param_name == "scaling_IgM")			Expt[i].Pop[j].Hosts[k].scaling_IgM = stod(param_value_string);
						if (param_name == "init_m")					Expt[i].Pop[j].Hosts[k].Params.init_m = stod(param_value_string);
						if (param_name == "init_g")					Expt[i].Pop[j].Hosts[k].Params.init_g = stod(param_value_string);
						if (param_name == "init_n")					Expt[i].Pop[j].Hosts[k].Params.init_n = stod(param_value_string);
						if (param_name == "psi")					Expt[i].Pop[j].Hosts[k].Params.psi = stod(param_value_string);
						if (param_name == "psi2")					Expt[i].Pop[j].Hosts[k].Params.psi2 = stod(param_value_string);
						if (param_name == "c")						Expt[i].Pop[j].Hosts[k].Params.c = stod(param_value_string);

	
						if (param_name == "hill_coeff")				Expt[i].Pop[j].Hosts[k].Params.hill_coeff = stod(param_value_string);
						if (param_name == "ec50")					Expt[i].Pop[j].Hosts[k].Params.ec50 = stod(param_value_string);
						if (param_name == "emax")					Expt[i].Pop[j].Hosts[k].PK.emax = stod(param_value_string);
						if (param_name == "cl_ed50")				Expt[i].Pop[j].Hosts[k].PK.cl_ed50 = stod(param_value_string);
						if (param_name == "cl_max")					Expt[i].Pop[j].Hosts[k].PK.cl_max = stod(param_value_string);
						if (param_name == "cl_var")					Expt[i].Pop[j].Hosts[k].PK.cl_var = stod(param_value_string);
						if (param_name == "cl_emax1")				Expt[i].Pop[j].Hosts[k].PK.cl_emax1 = stod(param_value_string);
						if (param_name == "cl_emax2")				Expt[i].Pop[j].Hosts[k].PK.cl_emax2 = stod(param_value_string);
						if (param_name == "V")						Expt[i].Pop[j].Hosts[k].PK.V = stod(param_value_string);
						if (param_name == "k_e")					Expt[i].Pop[j].Hosts[k].PK.k_e = stod(param_value_string);
						if (param_name == "R0")						Expt[i].Pop[j].Hosts[k].Params.R0= stod(param_value_string);
						
						if (param_name == "pk_model_number")		model.pk_model_number = stoi(param_value_string);
						if (param_name == "drug_moa")				model.drug_moa = stoi(param_value_string);
						if (param_name == "hierarchical")			model.hierarchical = stoi(param_value_string);
						if (param_name == "lag_g")					Expt[i].Pop[j].Hosts[k].Params.lag_g = stod(param_value_string);
						
						if (param_name == "no_fit")					Fit.no_fit = stoi(param_value_string);
						if (param_name == "no_particles")			Fit.no_particles = stoi(param_value_string);
						if (param_name == "iterations")				Fit.iterations = stod(param_value_string);
						if (param_name == "burn_in")				Fit.burn_in = stod(param_value_string);
						if (param_name == "thin")					Fit.thin = stoi(param_value_string);
						if (param_name == "scaling")				Fit.scaling = stod(param_value_string);
						if (param_name == "acpt_rate")				Fit.req_acpt_rate = stod(param_value_string);
						if (param_name == "no_post_samples")		Fit.no_post_samples = stoi(param_value_string);
						if (param_name == "cov_start")				Fit.cov_start = stoi(param_value_string);
						if (param_name == "cov_stop")				Fit.cov_stop = stoi(param_value_string);
						if (param_name == "cov_update")				Fit.cov_update = stoi(param_value_string);

						if (param_name == "error_v")				Expt[i].Pop[j].Hosts[k].Data.error_v = stod(param_value_string);
						if (param_name == "error_m")				Expt[i].Pop[j].Hosts[k].Data.error_m = stod(param_value_string);	
						if (param_name == "error_g")				Expt[i].Pop[j].Hosts[k].Data.error_g = stod(param_value_string);
						if (param_name == "error_n")				Expt[i].Pop[j].Hosts[k].Data.error_n = stod(param_value_string);
						if (param_name == "threshold")				model.threshold = stod(param_value_string);
						
				}
				
				Param_File.close();

			}
		}
	}
}

void Read_Hyperpriors(string filename, hyperprior& Hyperpriors, model_type model, expt* Expt) {

	gsl_rng_set(g2, long(time(NULL)));

	ifstream Hyper_priors;

	Hyper_priors.open(filename);

	string line;

	int no_par;

	no_par = 0;

	while (getline(Hyper_priors, line)) { ++no_par; }

	Hyper_priors.close();

	Hyperpriors.no_params = no_par;

	Hyperpriors.logscale = new int[no_par];

	Hyper_priors.open(filename);

	string param_name, distribution, par1, par2, logscale;

	int n = 0;

	while (!Hyper_priors.eof())
	{
	
		getline(Hyper_priors, param_name, '\t');
		getline(Hyper_priors, logscale, '\t');
		getline(Hyper_priors, distribution, '\n');
	
		Hyperpriors.Name.push_back(param_name);
		Hyperpriors.Dist.push_back(distribution);

		Hyperpriors.logscale[n] = stoi(logscale);
		
		n = n + 1;
	}

	Hyper_priors.close();

}

void Read_Param_Fit(string filename, int n_expt, fit_params& Fit, vector<string>& Param_Fit, vector<string>& Prior_Dist) {

	Fit.logscale = new int[Fit.no_fit];/*0-fitting parameter on linear scale, 1-fitting parameter on log10 scale*/
	Fit.level = new int[Fit.no_fit]; /*0-fitting at global level, 1- fitting at experiment level, 2-fitting at group level, 3- fitting at individual level*/
	Fit.block = new int[Fit.no_fit];
	Fit.init_sd = new double[Fit.no_fit];/*initial jump size for proposing values*/
	Fit.max_sd = new double[Fit.no_fit];/*max jump size for proposing values*/
	Fit.min_sd = new double[Fit.no_fit];/*min jump size for proposing values*/
	Fit.start_values = new double[Fit.no_fit];/*starting value of each param we want to fit*/
	Fit.sample = new int[Fit.no_fit];

	ifstream Par_Fit;

	Par_Fit.open(filename);

	string hyperparameter,logscale, init_sd, sd_max, sd_min, level, block, param_name, param_value_string, distribution, lower, upper, sample;

	int i =0;
	int n = 0;

	vector<string> temp3;
	vector<string> temp4;


	while (!Par_Fit.eof())
	{
		
		if (n == Fit.no_fit) break;

		getline(Par_Fit, param_name, '\t');
		getline(Par_Fit, param_value_string, '\t');
		getline(Par_Fit, hyperparameter, '\t');
		getline(Par_Fit, logscale, '\t');
		getline(Par_Fit, init_sd, '\t');
		getline(Par_Fit, sd_min, '\t');
		getline(Par_Fit, sd_max, '\t');
		getline(Par_Fit, level, '\t');
		getline(Par_Fit, sample, '\t');
		getline(Par_Fit, block, '\t');
		getline(Par_Fit, distribution, '\t');
		getline(Par_Fit, lower, '\t');
		getline(Par_Fit, upper, '\n');

		if (i == 0) {

			Param_Fit.push_back(param_name);
			Prior_Dist.push_back(distribution);
		}

		Fit.start_values[n] = stod(param_value_string);
		Fit.logscale[n] = stoi(logscale);
		Fit.init_sd[n] = stod(init_sd);
		Fit.max_sd[n] = stod(sd_max);
		Fit.min_sd[n] = stod(sd_min);
		Fit.level[n] = stoi(level);
		Fit.block[n] = stoi(block);
		Fit.sample[n] = stoi(sample);

		if (stoi(hyperparameter) == 0) {

			vector<string> temp1;
			vector<string> temp2;

			for (int j = 0; j < n_expt;j++) {

				temp1.push_back(lower);
				temp2.push_back(upper);
				
			}

			Fit.lower.push_back(temp1);
			Fit.upper.push_back(temp2);

			n = n + 1;
		}

		else {

			if (Param_Fit[n] == param_name) {i = i + 1;}

			temp3.push_back(lower);
			temp4.push_back(upper);

			if (i == n_expt) { 
				
				Fit.lower.push_back(temp3);
				Fit.upper.push_back(temp4);
				
				temp3.clear();
				temp4.clear();

				i = 0;
				n = n + 1; 			
				
			}

		}

	}

	Par_Fit.close();
}

void Read_Starting_Values(string filename, expt* Expt, model_type model, fit_params Fit) {

	ifstream Start_Values;

	string param_value;

	Start_Values.open(filename);

	for (int n = 0; n < Fit.no_fit; n++) {

		while (!Start_Values.eof()) {

			if (Fit.level[n] == 0) {

				getline(Start_Values, param_value, '\n');

				cout << param_value << endl;

				for (int i = 0; i < model.n_expt; i++) {

					for (int j = 0; j < Expt[i].no_groups; j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = stod(param_value);

							//cout << n<<','<<Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] << endl;
						}
					}
				}

			}

			if (Fit.level[n] == 1) {

				for (int i = 0; i < model.n_expt; i++) {

					getline(Start_Values, param_value, '\n');
					//cout << param_value << endl;

					for (int j = 0; j < Expt[i].no_groups; j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = stod(param_value);

							//cout << n << ',' << Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] << endl;
						}
					}
				}

			}

			if (Fit.level[n] == 2) {

				for (int i = 0; i < model.n_expt; i++) {

					for (int j = 0; j < Expt[i].no_groups; j++) {

						getline(Start_Values, param_value, '\n');
						//cout << param_value << endl;

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = stod(param_value);

							//cout << n << ',' << Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] << endl;
						}
					}
				}

			}

			if (Fit.level[n] == 3) {

				for (int i = 0; i < model.n_expt; i++) {

					for (int j = 0; j < Expt[i].no_groups; j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							getline(Start_Values, param_value, '\n');

							//cout << param_value << endl;

							Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = stod(param_value);

							//cout << n << ',' << Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] << endl;
						}
					}
				}

			}

			n = n + 1;
		}


		//system("pause");
	}

	Start_Values.close();
}

void Read_Dose_Schedule(string filename, expt& Expt) {

	ifstream Dose_Scheulde;
	Dose_Scheulde.open(filename);

	string hour, number;

	int no_doses;

	no_doses = -1;

	while (!Dose_Scheulde.eof()) {

		getline(Dose_Scheulde, number, '\t');
		getline(Dose_Scheulde, hour, '\n');

		no_doses += 1;

		if (Dose_Scheulde.eof()) break;
	}

	Dose_Scheulde.close();

	for (int j = 0; j < Expt.no_groups;j++) {

		for (int k = 0; k < Expt.Pop[j].no_hosts; k++) {

			Expt.Pop[j].Hosts[k].no_doses = no_doses;

			Expt.Pop[j].Hosts[k].Doses = new dose[no_doses];

		}
	}

	int t = 0;

	Dose_Scheulde.open(filename);

	while (!Dose_Scheulde.eof()) {

		getline(Dose_Scheulde, number, '\t');
		getline(Dose_Scheulde, hour, '\n');

		if (Dose_Scheulde.eof()) break;

		for (int j = 0; j < Expt.no_groups;j++) {

			for (int k = 0; k < Expt.Pop[j].no_hosts; k++) {

				Expt.Pop[j].Hosts[k].Doses[t].dose_no = stoi(number);
				Expt.Pop[j].Hosts[k].Doses[t].dose_time = stoi(hour);
				Expt.Pop[j].Hosts[k].Doses[t].dose_size = Expt.Pop[j].Hosts[k].dose_size;


				if (Expt.Pop[j].Hosts[k].Doses[t].dose_no == 0) { Expt.Pop[j].Hosts[k].first_dose = stoi(hour); }
				if (Expt.Pop[j].Hosts[k].Doses[t].dose_no == (Expt.Pop[j].Hosts[k].no_doses-1)) { Expt.Pop[j].Hosts[k].last_dose = stoi(hour); }


			}

		}

		t = t + 1;
	}

}

void Allocate_Memory(fit_params Fit, int n_expt, expt*Expt) {

	long length = long((Fit.iterations) / Fit.thin);

	long length_sub = long(((Fit.iterations - Fit.burn_in) / Fit.thin));

	for (int i = 0; i < n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int n = 0; n < Expt[i].Pop[j].no_hosts;n++) {

				Expt[i].Pop[j].Hosts[n].X = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].X[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) {
						Expt[i].Pop[j].Hosts[n].X[m][t] = new long double[RUNS];
					}
				}

				Expt[i].Pop[j].Hosts[n].Y = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].Y[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].Y[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].V = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].V[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].V[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].M = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].M[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].M[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].G = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles; m++) {
					Expt[i].Pop[j].Hosts[n].G[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS; t++) { Expt[i].Pop[j].Hosts[n].G[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].D = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].D[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].D[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].N = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].N[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].N[m][t] = new long double[RUNS]; }
				}
	

				Expt[i].Pop[j].Hosts[n].LL = new double* [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) { Expt[i].Pop[j].Hosts[n].LL[m] = new double[DAYS]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.C = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.C[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_C = new double* [length_sub];
				for (int k = 0; k < length_sub;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_C[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.LC = new double[length];
				Expt[i].Pop[j].Hosts[n].Par_Chain.R0 = new double[length];

				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_LC = new double[length_sub];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_R0 = new double[length_sub];

				Expt[i].Pop[j].Hosts[n].Par_Chain.SD = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.SD[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Posterior = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Posterior[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Prop_Acpt = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Prop_Acpt[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Param_Sample = new double* [Fit.no_post_samples];
				for (int k = 0; k < Fit.no_post_samples;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Param_Sample[k] = new double[Fit.no_fit]; }
			
				Expt[i].Pop[j].Hosts[n].Par_Chain.Mean_Param_Sample = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Par_Chain.W = new double* [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles; m++) { Expt[i].Pop[j].Hosts[n].Par_Chain.W[m] = new double[DAYS]; }
				
				Expt[i].Pop[j].Hosts[n].Par_Chain.Lower = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Median = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Upper = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_val = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_val = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_sd = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_sd = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_posterior = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_posterior = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_ll = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_ll = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_pr = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_pr = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.no_acpt = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Simulation = new host_predict[Fit.no_post_samples];
				for (int s = 0; s < Fit.no_post_samples;s++) {
					
					Expt[i].Pop[j].Hosts[n].Simulation[s].X = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].Y = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].V = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].M = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].G = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].D = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].N = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].R = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].Eff = new double[STEPS];
					
				}

				Expt[i].Pop[j].Hosts[n].Bounds = new host_predict[3];
				for (int s = 0; s < 3;s++) {

					Expt[i].Pop[j].Hosts[n].Bounds[s].X = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].Y = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].V = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].M = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].G = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].D = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].N = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].R= new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].Eff = new double[STEPS];
				
			
				}

				Expt[i].Pop[j].Hosts_Avg[n].X = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].X[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].Y = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].Y[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].V = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].V[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].M = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].M[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].G = new long double* [STEPS];
				for (int t = 0; t < STEPS; t++) {
					Expt[i].Pop[j].Hosts_Avg[n].G[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].D = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].D[t] = new long double[2];
				}

				Expt[i].Pop[j].Hosts_Avg[n].N = new long double* [STEPS];
				for (int t = 0;t < STEPS;t++) {
					Expt[i].Pop[j].Hosts_Avg[n].N[t] = new long double[2];
				}


				Expt[i].Pop[j].Hosts[n].PK.A1 = new double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].PK.A1[m] = new double* [STEPS];
					for (int t = 0; t < STEPS;t++) {
						Expt[i].Pop[j].Hosts[n].PK.A1[m][t] = new double[RUNS];
					}
				}

			}

			Expt[i].Pop[j].Avg.X = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.X[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.Y = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.Y[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.V = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.V[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.M = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.M[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.G = new long double* [STEPS];
			for (int t = 0; t < STEPS; t++) {
				Expt[i].Pop[j].Avg.G[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.D = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.D[t] = new long double[2];
			}

			Expt[i].Pop[j].Avg.N = new long double* [STEPS];
			for (int t = 0;t < STEPS;t++) {
				Expt[i].Pop[j].Avg.N[t] = new long double[2];
			}

		}
	}

	/*initialise arrays*/

	for (int i = 0; i < n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int n = 0; n < Expt[i].Pop[j].no_hosts; n++) {

				for (int t = 0; t < STEPS;t++) {

					for (int s = 0; s < 2;s++) {

						Expt[i].Pop[j].Hosts_Avg[n].X[t][s] = 0;
						Expt[i].Pop[j].Hosts_Avg[n].Y[t][s] = 0;
						Expt[i].Pop[j].Hosts_Avg[n].V[t][s] = 0;
						Expt[i].Pop[j].Hosts_Avg[n].M[t][s] = 0;
						Expt[i].Pop[j].Hosts_Avg[n].G[t][s] = 0;
						Expt[i].Pop[j].Hosts_Avg[n].D[t][s] = 0;
						
						Expt[i].Pop[j].Avg.X[t][s] = 0;
						Expt[i].Pop[j].Avg.Y[t][s] = 0;
						Expt[i].Pop[j].Avg.V[t][s] = 0;
						Expt[i].Pop[j].Avg.M[t][s] = 0;
						Expt[i].Pop[j].Avg.G[t][s] = 0; 
						Expt[i].Pop[j].Avg.D[t][s] = 0;
						Expt[i].Pop[j].Avg.N[t][s] = 0;


					}
				}
			}
		}
	}
}

void Init_Values(host& Host, int M, model_type model) {

	for (int m = 0; m < M;m++) {

		for (int r = 0; r < RUNS;r++) {

			/*initial values for each particle*/

			double diff = (1 - exp(-Host.Params.gamma*DT)) / (Host.Params.gamma*DT);

			Host.X[m][0][0] = Host.Params.A / (diff * Host.Params.gamma);

			Host.Y[m][0][0] = 0;
			Host.M[m][0][0] = Host.Params.init_m;
			Host.G[m][0][0] = Host.Params.init_g; 
			Host.D[m][0][0] = 0;
			Host.N[m][0][0] = Host.Params.init_n;

			Host.PK.A1[m][0][0] = 0;

			if (Host.time_infection==0 && RUN_IN==0) {

				Host.V[m][0][0] = Host.scaling_inoculum*Host.viral_inoculum;

			}

			else {

				Host.V[m][0][0] = 0;
			}

		}
	}
}

void Drug_Dosage(host &Host, model_type model) {

	int t1, D;

	for (int t = 0; t < STEPS;t++) { Host.dosage[t] = 0; }

	D = Host.no_doses;/*total number of doses*/

	for (int d = 0; d < D;d++) {

		t1 = Host.Doses[d].dose_time+24*RUN_IN;

		Host.dosage[t1] += Host.Doses[d].dose_size;

		}
}

double Cond_Dist(hyperprior& Hyperpriors, const vector<string>& Param_Fit, fit_params Fit, host& Host) {

	double val, a, b,result;
	
	int length, index;
	
	length = Hyperpriors.no_params;

	result = 0;

	for (int n = 0; n < length; n++) {

		if (strcmp(Hyperpriors.Name[n].c_str(), "beta") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "beta") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "beta"));

				if (Fit.logscale[index] == 0) { val = Host.Params.beta; }

				else { val = log10(Host.Params.beta); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.beta; }

				else { val = log10(Host.Params.beta); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "beta_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "beta_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.beta_par1; }

				else { a = log10(Host.Params.beta_par1); }

			}

			else {a = Host.Params.beta_par1;}

			if (find(Param_Fit.begin(), Param_Fit.end(), "beta_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "beta_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.beta_par2; }

				else { b = log10(Host.Params.beta_par2); }

			}

			else { b = Host.Params.beta_par2; }


			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2*gsl_ran_cauchy_pdf(val,b); }
			
		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "eta") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta"));

				if (Fit.logscale[index] == 0) { val = Host.Params.eta; }

				else { val = log10(Host.Params.eta); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.eta; }

				else { val = log10(Host.Params.eta); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.eta_par1; }

				else { a = log10(Host.Params.eta_par1); }

			}

			else { a = Host.Params.eta_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.eta_par2; }

				else { b = log10(Host.Params.eta_par2); }

			}

			else { b = Host.Params.eta_par2; }


			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "eta2") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta2"));

				if (Fit.logscale[index] == 0) { val = Host.Params.eta2; }

				else { val = log10(Host.Params.eta2); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.eta2; }

				else { val = log10(Host.Params.eta2); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta2_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta2_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.eta2_par1; }

				else { a = log10(Host.Params.eta2_par1); }

			}

			else {a = Host.Params.eta2_par1;}

			if (find(Param_Fit.begin(), Param_Fit.end(), "eta2_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "eta2_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.eta2_par2; }

				else { b = log10(Host.Params.eta2_par2); }

			}

			else { b = Host.Params.eta2_par2; }


			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "scaling_inoculum") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "scaling_inoculum") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "scaling_inoculum"));

				if (Fit.logscale[index] == 0) { val = Host.scaling_inoculum; }

				else { val = log10(Host.scaling_inoculum); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.scaling_inoculum; }

				else { val = log10(Host.scaling_inoculum); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "scaling_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "scaling_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.scaling_par1; }

				else { a = log10(Host.Params.scaling_par1); }

			}

			else {a = Host.Params.scaling_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "scaling_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "scaling_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.scaling_par2; }

				else { b = log10(Host.Params.scaling_par2); }

			}

			else {b = Host.Params.scaling_par2;}

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "scaling_IgM") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "scaling_IgM") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "scaling_IgM"));

				if (Fit.logscale[index] == 0) { val = Host.scaling_IgM; }

				else { val = log10(Host.scaling_IgM); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.scaling_IgM; }

				else { val = log10(Host.scaling_IgM); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "SFM_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "SFM_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.SFM_par1; }

				else { a = log10(Host.Params.SFM_par1); }

			}

			else {a = Host.Params.SFM_par1;}

			if (find(Param_Fit.begin(), Param_Fit.end(), "SFM_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "SFM_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.SFM_par2; }

				else { b = log10(Host.Params.SFM_par2); }

			}

			else {b = Host.Params.SFM_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }
	
		}
		
		if (strcmp(Hyperpriors.Name[n].c_str(), "psi") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi"));

				if (Fit.logscale[index] == 0) { val = Host.Params.psi; }

				else { val = log10(Host.Params.psi); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.psi; }

				else { val = log10(Host.Params.psi); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.psi_par1; }

				else { a = log10(Host.Params.psi_par1); }

			}

			else {a=Host.Params.psi_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.psi_par2; }

				else { b = log10(Host.Params.psi_par2); }

			}

			else {b = Host.Params.psi_par2;}

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }
	
		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "psi2") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi2"));

				if (Fit.logscale[index] == 0) { val = Host.Params.psi2; }

				else { val = log10(Host.Params.psi2); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.psi2; }

				else { val = log10(Host.Params.psi2); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi2_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi2_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.psi2_par1; }

				else { a = log10(Host.Params.psi2_par1); }

			}

			else {a = Host.Params.psi2_par1;}

			if (find(Param_Fit.begin(), Param_Fit.end(), "psi2_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "psi2_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.psi2_par2; }

				else { b = log10(Host.Params.psi2_par2); }

			}

			else {b = Host.Params.psi2_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "kappa_m") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "kappa_m") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "kappa_m"));

				if (Fit.logscale[index] == 0) { val = Host.Params.kappa_m; }

				else { val = log10(Host.Params.kappa_m); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.kappa_m; }

				else { val = log10(Host.Params.kappa_m); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "kappa_m_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "kappa_m_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.kappa_m_par1; }

				else { a = log10(Host.Params.kappa_m_par1); }

			}

			else { a = Host.Params.kappa_m_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "kappa_m_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "kappa_m_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.kappa_m_par2; }

				else { b = log10(Host.Params.kappa_m_par2); }

			}

			else { b = Host.Params.kappa_m_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "init_m") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_m") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_m"));

				if (Fit.logscale[index] == 0) { val = Host.Params.init_m; }

				else { val = log10(Host.Params.init_m); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.init_m; }

				else { val = log10(Host.Params.init_m); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_m_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_m_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.init_m_par1; }

				else { a = log10(Host.Params.init_m_par1); }

			}

			else { a = Host.Params.init_m_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_m_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_m_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.init_m_par2; }

				else { b = log10(Host.Params.init_m_par2); }

			}

			else { b = Host.Params.init_m_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "init_g") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_g") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_g"));

				if (Fit.logscale[index] == 0) { val = Host.Params.init_g; }

				else { val = log10(Host.Params.init_g); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.init_g; }

				else { val = log10(Host.Params.init_g); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_g_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_g_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.init_g_par1; }

				else { a = log10(Host.Params.init_g_par1); }

			}

			else { a = Host.Params.init_g_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "init_g_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "init_g_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.init_g_par2; }

				else { b = log10(Host.Params.init_g_par2); }

			}

			else { b = Host.Params.init_g_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}

		if (strcmp(Hyperpriors.Name[n].c_str(), "c") == 0) {

			if (find(Param_Fit.begin(), Param_Fit.end(), "c") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "cf"));

				if (Fit.logscale[index] == 0) { val = Host.Params.c; }

				else { val = log10(Host.Params.c); }

			}

			else {

				if (Hyperpriors.logscale[n] == 0) { val = Host.Params.c; }

				else { val = log10(Host.Params.c); }
			}

			if (find(Param_Fit.begin(), Param_Fit.end(), "c_par1") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "c_par1"));

				if (Fit.logscale[index] == 0) { a = Host.Params.c_par1; }

				else { a = log10(Host.Params.c_par1); }

			}

			else { a = Host.Params.c_par1; }

			if (find(Param_Fit.begin(), Param_Fit.end(), "c_par2") != Param_Fit.end()) {

				index = distance(Param_Fit.begin(), find(Param_Fit.begin(), Param_Fit.end(), "c_par2"));

				if (Fit.logscale[index] == 0) { b = Host.Params.c_par2; }

				else { b = log10(Host.Params.c_par2); }

			}

			else { b = Host.Params.c_par2; }

			if (strcmp(Hyperpriors.Dist[n].c_str(), "Uniform") == 0) { result += gsl_ran_flat_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Normal") == 0) { result += gsl_ran_gaussian_pdf(val - a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Gamma") == 0) { result += gsl_ran_gamma_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Beta") == 0) { result += gsl_ran_beta_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "Lognormal") == 0) { result += gsl_ran_lognormal_pdf(val, a, b); }
			if (strcmp(Hyperpriors.Dist[n].c_str(), "HalfCauchy") == 0) { result += 2 * gsl_ran_cauchy_pdf(val, b); }

		}
	}
	
	if (result == 0) { result = -50000; }

	else { result = log(result); }

	return result;
}

/*calculate mean across runs for each host, for each timestep*/
void Mean_Host(host_avg*Hosts_Avg, host*Host, int n, int runs, int steps) {

	for (int t = 0; t < steps;t++) {

		for (int r = 0; r < runs; r++) {

			Hosts_Avg[n].X[t][0] += Host[n].X[0][t][r] / double(runs);
			Hosts_Avg[n].Y[t][0] += Host[n].Y[0][t][r] / double(runs);
			Hosts_Avg[n].V[t][0] += Host[n].V[0][t][r] / double(runs);
			Hosts_Avg[n].N[t][0] += Host[n].N[0][t][r] / double(runs);
			Hosts_Avg[n].M[t][0] += Host[n].M[0][t][r] / double(runs);
			Hosts_Avg[n].G[t][0] += Host[n].G[0][t][r] / double(runs);
			Hosts_Avg[n].D[t][0] += Host[n].D[0][t][r] / double(runs);
		}
	}

}

/*calculate standard deviation between runs for each timestep*/
void Std_Host(host_avg*Hosts_Avg, host*Host, int n, int runs, int steps) {

	double	diff_x, diff_y, diff_v, diff_m,diff_g, diff_d, diff_n,
			var_x,var_y,var_v,var_m,var_g,var_d, var_n,
			sum_x,sum_y,sum_v,sum_m,sum_g, sum_d, sum_n;

	for (int t = 0; t < steps;t++) {

		sum_x = 0; sum_y = 0; sum_v = 0; sum_m = 0; sum_g = 0; sum_d = 0; sum_n = 0;

		for (int r = 0; r < runs; r++) {

			diff_x = Host[n].X[0][t][r] - Hosts_Avg[n].X[t][0];
			sum_x += pow(diff_x, 2);

			diff_y = Host[n].Y[0][t][r] - Hosts_Avg[n].Y[t][0];
			sum_y += pow(diff_y, 2);

			diff_v = Host[n].V[0][t][r] - Hosts_Avg[n].V[t][0];
			sum_v += pow(diff_v, 2);

			diff_m = Host[n].M[0][t][r] - Hosts_Avg[n].M[t][0];
			sum_m += pow(diff_m, 2);

			diff_g = Host[n].G[0][t][r] - Hosts_Avg[n].G[t][0];
			sum_g += pow(diff_g, 2);

			diff_d = Host[n].D[0][t][r] - Hosts_Avg[n].D[t][0];
			sum_d += pow(diff_d, 2);

			diff_n = Host[n].N[0][t][r] - Hosts_Avg[n].N[t][0];
			sum_n += pow(diff_n, 2);

		}

		var_x = sum_x / double(runs);
		Hosts_Avg[n].X[t][1] = sqrt(var_x);

		var_y = sum_y / double(runs);
		Hosts_Avg[n].Y[t][1] = sqrt(var_y);

		var_v = sum_v / double(runs);
		Hosts_Avg[n].V[t][1] = sqrt(var_v);

		var_m = sum_m / double(runs);
		Hosts_Avg[n].M[t][1] = sqrt(var_m);

		var_g = sum_g / double(runs);
		Hosts_Avg[n].G[t][1] = sqrt(var_g);

		var_d = sum_d / double(runs);
		Hosts_Avg[n].D[t][1] = sqrt(var_d);

		var_n = sum_n / double(runs);
		Hosts_Avg[n].N[t][1] = sqrt(var_n);
	}
}

void SubSample(double**C, double**S, double*LC, double*S2, double*R0, double*S3, int no_fit, long length, double burn_in, int thin) {

	long k2;

	int k0 = int(burn_in / thin);

	k2 = 0;

	for (int k = k0; k < length; k++) {

		S2[k2] = LC[k];
		S3[k2] = R0[k];

		for (int n = 0; n < no_fit; n++) {

			S[k2][n] = C[k][n];
		}

		k2 += 1;
	}
}

double Get_Quantile_Chains(double**C, int iterations, int n, double p) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	double*X;
	X = new double[iterations];

	for (int i = 0; i < iterations; i++) { X[i] = C[i][n]; }

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, p);

	return result;
}

double Get_Median(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.5);

	return result;
}

double Get_Lower(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.025);

	return result;
}

double Get_Upper(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.975);

	return result;
}

void Set_New_Parameter_Values(host &Host, const vector<string> &v, int no_fit, double*S) {

	for (int n = 0; n < no_fit; n++) {

		if (strcmp(v[n].c_str(), "beta") == 0) { Host.Params.beta = S[n]; }
		if (strcmp(v[n].c_str(), "kappa") == 0) { Host.Params.kappa = S[n]; }
		if (strcmp(v[n].c_str(), "kappa_m") == 0) { Host.Params.kappa_m = S[n]; }
		if (strcmp(v[n].c_str(), "kappa_n") == 0) { Host.Params.kappa_n = S[n]; }
		if (strcmp(v[n].c_str(), "rho") == 0) { Host.Params.rho = S[n]; }
		if (strcmp(v[n].c_str(), "eta") == 0) { Host.Params.eta = S[n]; }
		if (strcmp(v[n].c_str(), "eta2") == 0) { Host.Params.eta2 = S[n]; }
		if (strcmp(v[n].c_str(), "ec50") == 0) { Host.Params.ec50 = S[n]; }
		if (strcmp(v[n].c_str(), "hill_coeff") == 0) { Host.Params.hill_coeff = S[n]; }
		if (strcmp(v[n].c_str(), "scaling_inoculum") == 0) { Host.scaling_inoculum = S[n]; }
		if (strcmp(v[n].c_str(), "init_m") == 0) { Host.Params.init_m = S[n]; }
		if (strcmp(v[n].c_str(), "init_g") == 0) { Host.Params.init_g = S[n]; }
		if (strcmp(v[n].c_str(), "beta_par1") == 0) { Host.Params.beta_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "beta_par2") == 0) { Host.Params.beta_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "eta_par1") == 0) { Host.Params.eta_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "eta_par2") == 0) { Host.Params.eta_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "eta2_par1") == 0) { Host.Params.eta2_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "eta2_par2") == 0) { Host.Params.eta2_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "scaling_par1") == 0) { Host.Params.scaling_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "scaling_par2") == 0) { Host.Params.scaling_par2= S[n]; }
		if (strcmp(v[n].c_str(), "scaling_IgM") == 0) { Host.scaling_IgM = S[n]; }
		if (strcmp(v[n].c_str(), "SFM_par1") == 0) { Host.Params.SFM_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "SFM_par2") == 0) { Host.Params.SFM_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "psi") == 0) { Host.Params.psi = S[n]; }
		if (strcmp(v[n].c_str(), "psi2") == 0) { Host.Params.psi2 = S[n]; }
		if (strcmp(v[n].c_str(), "psi_par1") == 0) { Host.Params.psi_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "psi_par2") == 0) { Host.Params.psi_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "psi2_par1") == 0) { Host.Params.psi2_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "psi2_par2") == 0) { Host.Params.psi2_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "kappa_m_par1") == 0) { Host.Params.kappa_m_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "kappa_m_par2") == 0) { Host.Params.kappa_m_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "init_m_par1") == 0) { Host.Params.init_m_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "init_m_par2") == 0) { Host.Params.init_m_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "init_g_par1") == 0) { Host.Params.init_g_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "init_g_par2") == 0) { Host.Params.init_g_par2 = S[n]; }
		if (strcmp(v[n].c_str(), "lag_g") == 0) { Host.Params.lag_g = S[n]; }
		if (strcmp(v[n].c_str(), "error_v") == 0) { Host.Data.error_v = S[n]; }
		if (strcmp(v[n].c_str(), "error_m") == 0) { Host.Data.error_m = S[n]; }
		if (strcmp(v[n].c_str(), "error_n") == 0) { Host.Data.error_n = S[n]; }
		if (strcmp(v[n].c_str(), "error_g") == 0) { Host.Data.error_g = S[n]; }
		if (strcmp(v[n].c_str(), "c") == 0) { Host.Params.c = S[n]; }
		if (strcmp(v[n].c_str(), "c_par1") == 0) { Host.Params.c_par1 = S[n]; }
		if (strcmp(v[n].c_str(), "c_par2") == 0) { Host.Params.c_par2 = S[n]; }

	}
}

void Get_Quantile_Simulations(host& Host,host_predict*Simulation,host_predict Bound, int no_simulations, int T, double p){

	double temp;
	double*S;
	S = new double[no_simulations];

	for (int t = 0; t < T; t++) {

		/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].X[t] > Simulation[j].X[t]) {

						temp = Simulation[i].X[t];
						Simulation[i].X[t] = Simulation[j].X[t];
						Simulation[j].X[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].X[t]; }

		/*get quatile at time t*/
			
			Bound.X[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);
	
		/*sort values in ascending order*/
			
			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].Y[t] > Simulation[j].Y[t]) {

						temp = Simulation[i].Y[t];
						Simulation[i].Y[t] = Simulation[j].Y[t];
						Simulation[j].Y[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].Y[t]; }
			
		/*get quatile at time t*/
			
			Bound.Y[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);
			
		/*sort values in ascending order*/
			
			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].V[t] > Simulation[j].V[t]) {

						temp = Simulation[i].V[t];
						Simulation[i].V[t] = Simulation[j].V[t];
						Simulation[j].V[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].V[t];}
			
		/*get quatile at time t*/
			
			Bound.V[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/
			
			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].M[t] > Simulation[j].M[t]) {

						temp = Simulation[i].M[t];
						Simulation[i].M[t] = Simulation[j].M[t];
						Simulation[j].M[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].M[t]; }
			
		/*get quatile at time t*/
			
			Bound.M[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].G[t] > Simulation[j].G[t]) {

						temp = Simulation[i].G[t];
						Simulation[i].G[t] = Simulation[j].G[t];
						Simulation[j].G[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].G[t]; }

			/*get quatile at time t*/

			Bound.G[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/
			
			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].D[t] > Simulation[j].D[t]) {

						temp = Simulation[i].D[t];
						Simulation[i].D[t] = Simulation[j].D[t];
						Simulation[j].D[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].D[t]; }
		
		/*get quatile at time t*/
			
			Bound.D[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

			/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].N[t] > Simulation[j].N[t]) {

						temp = Simulation[i].N[t];
						Simulation[i].N[t] = Simulation[j].N[t];
						Simulation[j].N[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].N[t]; }

			/*get quatile at time t*/

			Bound.N[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);


			/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].R[t] > Simulation[j].R[t]) {

						temp = Simulation[i].R[t];
						Simulation[i].R[t] = Simulation[j].R[t];
						Simulation[j].R[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].R[t]; }

			/*get quatile at time t*/

			Bound.R[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

			/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].Eff[t] > Simulation[j].Eff[t]) {

						temp = Simulation[i].Eff[t];
						Simulation[i].Eff[t] = Simulation[j].Eff[t];
						Simulation[j].Eff[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].Eff[t]; }

			/*get quatile at time t*/

			Bound.Eff[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

	}
}

void Find_Mean_Value(double**X, double *Y, int n, int m) {

	for (int i = 0; i < n; i++) {

		Y[i] = 0;

		for (int j = 0; j < m; j++) {

			Y[i] += X[j][i]/m;
		}

	}

}

double Calculate_Rt(host& Host, int t, model_type model) {

	double Rt, effect, conc;
	int t1;

	t1 = t;

	conc = Host.D[0][t][0];

	effect = Host.PK.emax / (1 + pow((Host.Params.ec50 / conc), Host.Params.hill_coeff));

		double prod_rate;

		if (Host.V[0][t1][0] < model.threshold) { prod_rate = 0; }

		else { prod_rate = Host.Params.omega; }

		Rt = (Host.Params.beta * Host.X[0][t1][0] * prod_rate * (1 - effect)) / (Host.Params.kappa * (Host.Params.delta + Host.Params.alpha * (1/(1+Host.Params.c*Host.M[0][t1][0]))* Host.M[0][t1][0]));
	
	return Rt;
}

void Transform_Values(expt* Expt, fit_params Fit, int length, model_type model) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				for (int l = 0; l < length; l++) {

					for (int n = 0; n < Fit.no_fit; n++) {

						if (Fit.logscale[n] == 1) { Expt[i].Pop[j].Hosts[k].Par_Chain.C[l][n] = pow(10, Expt[i].Pop[j].Hosts[k].Par_Chain.C[l][n]); }

					}
				}
			}
		}
	}
}

void Calculate_Quantiles(expt* Expt, fit_params Fit, model_type model, int no_samples) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Median = Get_Median(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Median = Get_Median(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Lower = Get_Lower(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Lower = Get_Lower(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Upper = Get_Upper(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Upper = Get_Upper(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				for (int n = 0; n < Fit.no_fit; n++) {

					Expt[i].Pop[j].Hosts[k].Par_Chain.Median[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.5);
					Expt[i].Pop[j].Hosts[k].Par_Chain.Lower[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.025);
					Expt[i].Pop[j].Hosts[k].Par_Chain.Upper[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.975);

				}
			}
		}
	}

}

void Posterior_Sample(expt*Expt, fit_params Fit, int iterations, model_type model) {

	int x;
	double* W;
	W = new double[iterations];

	for (int k = 0; k < iterations; k++) {  W[k] = 1 / iterations; }

	/*all values equally likely to be sampled*/

	F1 = gsl_ran_discrete_preproc(iterations, W);
	
	for (int s = 0; s < Fit.no_post_samples; s++) {

		x = int(gsl_ran_discrete(g2, F1));

		for (int n = 0; n <Fit.no_fit; n++) {

			switch(Fit.level[n]){

			case 0:

				//x = int(gsl_ran_discrete(g2, F1));

			//	if (n == 0) { cout << s<<','<<x << endl; }

				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[0].Pop[0].Hosts[0].Par_Chain.Sub_C[x][n];
							

						}
					}
				}

				//if (n == 0) { cout << s << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.Param_Sample[s][n] << endl; }

				break;

			case 1:

				for (int i = 0; i < model.n_expt;i++) {

					//x = int(gsl_ran_discrete(g2, F1));

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[i].Pop[0].Hosts[0].Par_Chain.Sub_C[x][n];
						}
					}
				}

				break;

			case 2:
		
				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						//x = int(gsl_ran_discrete(g2, F1));

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[i].Pop[j].Hosts[0].Par_Chain.Sub_C[x][n];
						}
					}
				}

				break;

			case 3:

				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

						//	x = int(gsl_ran_discrete(g2, F1));

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C[x][n];
						}
					}
				}

				break;
			}
		}
	}

	//for (int s = 0; s < Fit.no_post_samples; s++) {

	//	for (int n = 0; n < Fit.no_fit; n++) {

	//		cout << s << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.Param_Sample[s][n] << endl;

	//	}
	//}
	//system("pause");

	gsl_ran_discrete_free(F1);

	delete[] W;
}

double Calculate_DIC(expt*Expt,fit_params Fit, model_type model, const vector<string>&Param_Fit, hyperprior& Hyperpriors) {

	double pd, dic;

	pd = 0; dic = 0;

	/*For each estimated parameter, find mean value of posterior samples*/

		for (int i = 0; i < model.n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

					Find_Mean_Value(Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample, Expt[i].Pop[j].Hosts[k].Par_Chain.Mean_Param_Sample, Fit.no_fit, Fit.no_post_samples);
				}
			}
		}


	/*Calculate deviance at mean posterior values*/

		double dev_post_mean, dev_post_sample, mean_dev_post_sample, ll;

		dev_post_mean = 0;
		ll = 0;

		for (int i = 0; i < model.n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

					Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Mean_Param_Sample);

					ll=Exact_Likelihood(Fit, 0, 0, Expt[i].Pop[j].Hosts[k], Param_Fit, Hyperpriors, model);

					dev_post_mean += -2 * ll;

				}
			}
		}

	/*Calculate mean of deviance evaluted at each individual posterior sample*/

		mean_dev_post_sample = 0;

		for (int s = 0; s < Fit.no_post_samples; s++) {

			dev_post_sample = 0;

			ll = 0;

			for (int i = 0; i < model.n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

						Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s]);

						ll=Exact_Likelihood(Fit, 0, 0, Expt[i].Pop[j].Hosts[k], Param_Fit, Hyperpriors, model);
			
						dev_post_sample += -2 * ll;
					}
				}
			}

			mean_dev_post_sample += dev_post_sample / Fit.no_post_samples;
		}

	/*Calculate DIC*/

		pd += mean_dev_post_sample - dev_post_mean;

		dic += mean_dev_post_sample + pd;

	return dic;
}

double Block_Size(expt* Expt, fit_params Fit, model_type model) {

	double size=0;

	for (int n = 0; n < Fit.no_fit; n++) {

		if (Fit.block[n] == 1) {

			if (Fit.level[n] == 0) { size += 1; }

			if (Fit.level[n] == 1) { size += 1*model.n_expt; }

			if (Fit.level[n] == 2) {

				for (int i = 0; i < model.n_expt;i++) {

					size += Expt[i].no_groups;
				}

			}

			if (Fit.level[n] == 3) {

				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						size += Expt[i].Pop[j].no_hosts;
					}
				}
			}
		}
	}

	return size;
}

void Covariance_Matrix(int size, int k1, int k2, expt* Expt, fit_params Fit, model_type model, gsl_matrix* Cov, gsl_matrix* Corr) {

	double corr, cov;

	double** temp; /*stores parameter values between iterations k1 and k2*/
	temp = new double* [size];

	int m = k2 - k1;

	for (int s = 0; s < size; s++) { temp[s] = new double[m]; }

	int s = 0;


	for (int n = 0; n < Fit.no_fit; n++) {

		switch (Fit.level[n]) {

		case 0:

			for (int k = k1; k < k2; k++) {
	
				temp[s][k - k1] = Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n];
				
				
			}

			s = s + 1;

			break;

		case 1:

			for (int i = 0; i < model.n_expt;i++) {

				for (int k = k1; k < k2; k++) {

					temp[s][k - k1] = Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n];
				}

				s = s + 1;
			}

			break;

		case 2:

			for (int i = 0; i < model.n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int k = k1; k < k2; k++) {

						temp[s][k - k1] = Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n];
					}

					s = s + 1;
				}
			}

			break;

		case 3:

			for (int i = 0; i < model.n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int l = 0; l < Expt[i].Pop[j].no_hosts;l++) {

						for (int k = k1; k < k2; k++) {

							temp[s][k - k1] = Expt[i].Pop[j].Hosts[l].Par_Chain.C[k][n];


						}

						s = s + 1;
					}
				}
			}

			break;

		}

	}


	/*compute covariance/correlation matrix*/

	for (int s1 = 0; s1 < size; s1++) {

		for (int s2 = 0; s2 < size; s2++) {

			gsl_matrix_set(Corr, s1, s2, 0); /*reset correlation matrix to zero*/
			corr = gsl_stats_correlation(temp[s1], 1, temp[s2], 1, m); /*computes correlation efficient between temp[n] and temp[m]*/
			gsl_matrix_set(Corr, s1, s2, corr); /*populate correlation matrix*/

			gsl_matrix_set(Cov, s1, s2, 0);/*reset covariance matrix to zero*/
			cov = gsl_stats_covariance(temp[s1], 1, temp[s2], 1, m);/*computes covariance between temp[n] and temp[m]*/
			gsl_matrix_set(Cov, s1, s2, cov);/*populate covariance matrix*/

		}
	}

	for (int s = 0; s < size; s++) { delete temp[s]; }

	delete[] temp;
}

