#include "matrix.h"
#include "data.h"
#include "solve.h"
#include <math.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES

#ifndef EPF_H
#define EPF_H

using namespace std;

/*---------------------------------------
			Code chunks
---------------------------------------*/


void sim(dynprob &prob, matrix &randus, data &dat);

imat randi(int n, int m);

matrix report_all_mc_out(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus);

matrix report_all_mc(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus);
matrix report_all_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);


matrix report_all_out_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);

imat permute(int p, int n, int no);

ord quickOrder(matrix arrr, imat orders, int left, int right);

matrix quickSort(matrix arr, int left, int right);

double cmedian(matrix vecs);

matrix var_infl(data dat, matrix (*infl) (data));
int find_nearest(double x, matrix vect);

void print_sim(data dat, string file);

void print_mat(matrix mat, string file);

void app_mat(matrix mat, string file, int seed);

void app_rmat(matrix mat, string file, int seed, matrix f);

void print_value(dynprob vsol, string prefx);

matrix readcsv(string files, int rows, int cols);

matrix mom_tm(data dat);

matrix mom_tm_cv(data dat);

double cv_epan(matrix &x, double h);

void band_epan(matrix &x, double fb, double a, double b, double c, double tol, double &ret);

void print_imat(imat mat, string file);

matrix mom_tm_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);

matrix mom_tm_cv_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);

matrix mom_tm_est_noint(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);


double normalCDF(double value);

double normalINV(double p);

matrix mom_tmINFL_nodemean(data dat);

void readdat(string files, int nobs, data &dat);

matrix epf_quad(data &dat, int degree);
matrix epf_quad_trans(data &dat, int degree);

matrix epf_cube(data &dat, int degree);

matrix epf_quad_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);
matrix epf_quad_trans_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);
matrix epf_cube_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);

matrix mom_tmINFL(data dat);
matrix mom_tm_cvINFL(data dat);

matrix epf_quadINFL(data dat);
matrix epf_quad_transINFL(data dat);

matrix epf_cubeINFL(data dat);

matrix test_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat);

#endif
