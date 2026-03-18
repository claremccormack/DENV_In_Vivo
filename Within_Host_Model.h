#ifndef WITHIN_HOST_MODEL_H
#define WITHIN_HOST_MODEL_H

#include "Structures.h"

void Within_Host_Model(host &Host, int thread, int start, int steps, int runs, model_type model);
long double IgM_Response(host Host, int thread, int t, int r, model_type model, double p1, double p2);
long double IgG_Response(host Host, int thread, int t, int r, model_type model, double p1, double p2);
long double Virus_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
double PK_Model(host &Host, int thread, int t, int r, model_type model);
double Drug_Absorption(host &Host, int t);
double Infection_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
long double Antigen_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
#endif
