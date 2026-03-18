#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Structures.h"
#include <vector>
#include <gsl/gsl_matrix.h>

using namespace std;

void Read_Data(string filename, expt*&Expt, int&n_expt);
void Read_Param_Values(string filename, model_type& model, fit_params& Fit, expt* Expt);
void Read_Hyperpriors(string filename, hyperprior &Hyperpriors, model_type model, expt* Expt);
void Read_Param_Fit(string filename, int n_expt, fit_params& Fit, vector<string>& Param_Fit, vector<string>& Param_Dist);
void Read_Dose_Schedule(string filename, expt&Expt);
double Cond_Dist(hyperprior& Hyperpriors, const vector<string>& Param_Fit, fit_params Fit, host& Host);
void Allocate_Memory(fit_params Fit, int n_expt, expt*Expt);
void Init_Values(host& Host, int M, model_type model);
void Drug_Dosage(host& Host, model_type model);
void Mean_Host(host_avg*Host_Avg, host*Host, int n, int runs, int steps);
void Std_Host(host_avg*Host_Avg, host*Host, int n, int runs, int steps);
void SubSample(double**C, double**S, double*LC, double*S2, double*R0, double*S3, int no_fit, long length, double burn_in, int thin);
double Get_Quantile_Chains(double**C, int iterations, int n, double p);
double Get_Median(double*C, int iterations);
double Get_Lower(double*C, int iterations);
double Get_Upper(double*C, int iterations);
void Set_New_Parameter_Values(host &Host, const vector<string> &v, int no_fit, double*S);
void Get_Quantile_Simulations(host& Host,host_predict*Simulation, host_predict Bound, int no_simulations, int T, double p);
void Find_Mean_Value(double**X, double*Y, int n, int m);
void Transform_Values(expt* Expt, fit_params Fit, int length, model_type model);
void Calculate_Quantiles(expt* Expt, fit_params Fit, model_type model, int no_samples);
void Posterior_Sample(expt* Expt, fit_params Fit, int iterations, model_type model);
double Calculate_DIC(expt* Expt, fit_params Fit, model_type model, const vector<string>& Param_Fit, hyperprior& Hyperpriors);
double Block_Size(expt* Expt, fit_params Fit, model_type model);
void Covariance_Matrix(int size, int k1, int k2, expt* Expt, fit_params Fit, model_type model, gsl_matrix* Cov, gsl_matrix* Corr);
double Calculate_Rt(host& Host, int t, model_type model);
void Read_Starting_Values(string filename, expt* Expt, model_type model, fit_params Fit);
#endif


