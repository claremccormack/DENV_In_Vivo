#ifndef STRUCTURES_H
#define STRUCTURES_H

#define RUNS 1
#define RUN_IN 0
#define DAYS 29+RUN_IN
#define DT .04166667 /*timestep*/
#define STEPS 24*(DAYS)+2

#include <vector>
#include <string>

using namespace std;

typedef struct MODEL_TYPE {

	int hierarchical; /*0-non-hierarchical, 1-hierarchical*/
	int pk_model_number; /*allow for different pharmacokinetic model structures*/
	int drug_moa; /*allow for different mechanisms of action*/
	
	double threshold;/*threshold for virion production to occur*/

	//int clearance_function_type;
	int run_mcmc;
	int n_expt; /*number of experiments*/
	std::string param_adjust;

}model_type;

typedef struct MODEL_PARAMS {

	/*params for in-vivo model*/

		double A; /*target cell production rate*/
		double alpha; /*infected cell clearance rate*/
		double beta; /*target cell infection rate (per virion)*/
		double delta; /*infected cell mortality rate*/
		double gamma; /*target cell mortality rate*/
		double kappa; /*virus clearance rate*/
		double omega; /*virus proliferation rate*/
		double c;

		double init_m; /*initial number of IgM*/
		double init_g; /*initial number of IgM*/
		double init_n; /*initial number of NS1*/
		double R0; /*basic reproduction number*/

		double eta; /*IgM response rate*/
		double eta2;/*IgM threshold parameter*/
		double psi;/*IgG response rate*/
		double psi2;/*IgG threshold parameter*/
		double rho; /*NS1 proliferation rate*/

		double kappa_n; /*NS1 decay rate*/
		double kappa_m; /*IgM decay rate*/
		double kappa_g; /*IgG decay rate*/

		double lag_g;

		double hill_coeff;
		double ec50;

	/*hyperparameters*/

		double beta_par1;
		double beta_par2;

		double eta_par1;
		double eta_par2;

		double eta2_par1;
		double eta2_par2;

		double scaling_par1;
		double scaling_par2;
		
		double kappa_m_par1;
		double kappa_m_par2;

		double init_m_par1;
		double init_m_par2;

		double init_g_par1;
		double init_g_par2;
		
		double psi_par1;
		double psi_par2;

		double psi2_par1;
		double psi2_par2;

		double SFM_par1;
		double SFM_par2;

		double c_par1;
		double c_par2;

}model_params;

typedef struct HYPERPRIOR {

	std::vector<std::string>Name;
	std::vector<std::string>Dist;
	int*logscale;
	int no_params;

}hyperprior;

typedef struct DOSE {

	int dose_no;
	int dose_time;
	double dose_size;

}dose;

typedef struct PK {

	double*** A1; /*first compartment*/

	/*parameters for Janseen clearance model*/

	double cl_ed50; /*dose which achieves 50% of max clearance*/
	double cl_emax1; /*max clearance for one dose per day*/
	double cl_emax2; /*max clearance for two doses per day*/
	double cl_max; /*maximum possible clearance*/
	double cl_var; /*variation in clearance between individuals*/
	double auc; /*area under the curve*/
	double k_e;/*elimination rate*/

	/*parameters for concentration-effect model*/

	double emax;
	double V;/*volume of distribution*/
	double k_a; /*absorption rate*/



}pk;

typedef struct HOST_DATA {

	int no_measurement_times;
	int Time[STEPS];

	/*limits of detection and measurement error levels*/

		double lod_v; /*limit of detection (virus)*/
		double lod_m; /*limit of detection (IgM)*/
		double lod_g; /*limit of detection (IgG)*/
		double lod_n; /*limit of detection (NS1)*/

		double error_v; /*measurement error (virus)*/
		double error_m; /*measurement error (IgM)*/
		double error_g; /*measurement error (IgG)*/
		double error_n; /*measurement error (NS1)*/

	/*arrays to store data*/

		double V[STEPS]; /*extracellular virus*/
		double D[STEPS]; /*drug concentration*/
		double N[STEPS]; /*NS1 concentration*/
		double M[STEPS]; /*IgM antibodies*/
		double G[STEPS]; /*IgG antibodies*/

}host_data;

typedef struct HOST_AVG {

	/*arrays to stroe average values*/

		long double**X; /*target cells*/
		long double**Y; /*infectious cells*/
		long double**V; /*virus*/
		long double**D; /*drug*/
		long double**N; /*NS1 antigen*/
		long double**M; /*IgM*/
		long double**G;/*IgG*/

}host_avg;

typedef struct HOST_PREDICT {

	/*stores predicted values*/

	long double* X; /*target cells*/
	long double* Y; /*infectious cells*/
	long double* V; /*virus*/
	long double* D; /*drug*/
	long double* N; /*ns1 antigen*/
	long double* M; /*IgM antigen*/
	long double* G; /*IgG antigen*/
	long double* R;
	double* Eff;

}host_predict;

typedef struct PARAM_CHAIN {

	double*curr_val; /*current parameter value*/
	double*new_val; /*proposed parameter value*/
	double*curr_ll; /*log-likelihood value based of current parameter values*/
	double*new_ll; /*log-likelihood value based of proposed parameter values*/
	double*curr_pr; /*prior value based of current parameter values*/
	double*new_pr; /*prior value based of proposed parameter values*/
	double*curr_posterior; /*posterior value based of current parameter values*/
	double*new_posterior; /*posterior value based on proposed parameter values*/
	double*curr_sd; /*current jump size*/
	double*new_sd; /*new jump size*/
	double*no_acpt; /*number of accepted values*/
	double curr_R0; /*R0 value based on current parameter values*/
	double new_R0; /*R0 value based on proposed parameter values*/

	double**C; /*chain of estimated parameter values*/
	double*LC; /*chain of log-likelihood values*/
	double*R0; /*chain of R0 values*/
	double**Sub_C; /*subset of estimated parameter values*/
	double*Sub_LC; /*subset of log-likelihood value*/
	double*Sub_R0; /*subset of R0 values*/
	double**SD; /*chain of jump sizes*/
	double**Posterior; /*posterior values*/
	double**Prop_Acpt; /*proportion of proposed parameter values accepted*/

	/*credible intervals*/

		double*Lower; /*lower bound of 95% credible interval for estimated parameter values*/
		double*Median; /*median estimated parameter values*/
		double*Upper; /*upper bound of 95% credible interval for estimated parameter values*/
		double LL_Median; /*median estimated log-likelihood value*/
		double R0_Median; /*median estimated R0 value*/
		double LL_Lower; /*lower bound of 95% credible interval for log-likelihood values*/
		double R0_Lower; /*lower bound of 95% credible interval for R0 values*/
		double LL_Upper; /*upper bound of 95% credible interval for log-likelihood values*/
		double R0_Upper; /*upper bound of 95% credible interval for R0 values*/

	double**W;
	double**Param_Sample;
	double*Mean_Param_Sample;
	
}param_chain;

typedef struct FIT_PARAMS {

	int no_fit; /*number of parameters we want to estimate*/
	double iterations; /*number of MCMC iterations*/
	double burn_in; /*length of burn-in for each MCMC chain*/
	int no_post_samples; /*number of samples from posterior parameter distributions*/
	double req_acpt_rate; /*desired MCMC accpetance rate*/
	int block_update; /*0-sequentially update parameter values, 1-block update parameters*/
	int thin; /*factor by which we thin MCMC chains*/

	double*start_values; /*initial value of each parameter we're estimating*/
	int*logscale; /*0- fit parameter on linear scale, 1- fit parameter on log10 scale*/
	int*level; /*0-fitting parameter globally, 1- fitting parameter locally to each experiment, 2-fitting parameter locally to each group, within an experiment, 3- fitting parameter to each individual host;*/
	int* block;
	vector<vector<string>>lower; /*lower bound for prior distribution*/
	vector<vector<string>>upper; /*upper bound for prior distribution*/
	int*sample;
	double*init_sd; /*initial jump size for each parameter we're estimating*/
	double*max_sd; /*max jump size for each parameter we're estimating*/
	double*min_sd; /*min jump size for each parameter we're estimating*/
	double scaling; /*controls how quickly we adjust jump size - large values give smaller, slower adjustments*/

	int no_particles;
	int trajectory;
	int cov_start;
	int cov_stop;
	int cov_update;

}fit_params;

/*hierarchical structure: experiment-group-individual host*/

typedef struct HOST {

	std::string id;
	model_params Params; /*viral kinetics model parameters*/
	pk PK; /*pk model parameters*/
	host_data Data; /*data*/

	param_chain Par_Chain; /*results of MCMC*/
	host_predict*Simulation; /*example simulations using parameter values estimated during MCMC*/
	host_predict*Bounds; /*stores posterior predictive intervals*/

	/*basic info about each host*/
	
		double viral_inoculum;
		double scaling_inoculum;
		double scaling_IgM;
		int time_infection;

		int no_doses;
		int first_dose;
		int last_dose;
		double dose_size;
		double dose_frequency;
		double dosage[STEPS];	
		dose* Doses;
		
	/*arrays to store cell populations for each host*/

		long double***X; /*target cells*/
		long double***Y; /*infectious cells*/
		long double***V; /*virus*/
		long double***N; /*NS1 antigen*/
		long double***D; /*drug*/
		long double***M; /*IgM*/
		long double***G; /*IgG*/

	double** LL;

}host;

typedef struct POP {

	int no_hosts; /*number of hosts within a group (within an experiment)*/
	host*Hosts;
	host_avg*Hosts_Avg;
	host_avg Avg;

}pop;

typedef struct EXPT {

	int no_groups; /*number of groups within an experiment*/
	pop* Pop;
	
}expt;

#endif
