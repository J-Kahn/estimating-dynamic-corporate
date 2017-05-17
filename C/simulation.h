#include "matrix.h"
#include "solve_tax.h"
#include "data.h"

#ifndef SIM_H
#define SIM_H
double normalINV( double p );

matrix qmarkov(matrix randu, double init, matrix A, matrix P);

matrix qdist(matrix z, double prev, double phi);

data sim(dynprob &prob, matrix &randus, matrix &randur, const int nbr);

data sim_natural(dynprob &prob, matrix &randus, matrix &randur, const int nbr);

data sim_markov(dynprob &prob, matrix &randus, matrix &randur, const int nbr);

data sim_permanent(dynprob &prob, matrix &randus, matrix &randur, const int nbr, const int nchange);

data sim_constant(dynprob &prob, matrix &randus, matrix &randur, const int nbr, int turn);

matrix moments(data dat);

matrix moments_natural(data dat);


matrix moments_est(dynprob &prob, data &ret, matrix &randus, matrix &randur);

data policy_data(dynprob &prob);

matrix compstats(data dat, matrix &ret1, matrix &ret2);

matrix compstats_natural(data dat, matrix &ret1, matrix &ret2, matrix &ret3, matrix &ret4, matrix &ret5, matrix &ret6);

matrix compstat_est(dynprob &prob, data &ret, matrix &randus, matrix &randur, matrix &ret1, matrix &ret2, matrix &ret3, matrix &ret4, matrix &ret5, matrix &ret6);
#endif