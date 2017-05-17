#include <math.h>
#include <iostream>
#include <fstream>
#include "epf.h"
#include "solve.h"
#define _USE_MATH_DEFINES

#ifndef OPT_H
#define OPT_H


matrix report_all_mc(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus);

void parameter_update(matrix &par, dynprob &prob);
void parameter_update(dynprob &prob);


double objec_mom(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_mom_cv(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_mom_mod(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_mom_kt(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_mom_hp(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

matrix NelderMead(matrix init_params, int maxevals, double tol, matrix m, matrix W, int seedn, matrix &f,
                      matrix minb, matrix maxb, dynprob &vali, matrix &randus, int nob, string des,
                      double (*objec) (matrix, matrix, matrix, dynprob&, data&, matrix&, int, matrix&));

double objec_epf_bspline(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_epf_quad(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);
double objec_epf_cube(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_epf_quad_trans(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);

double objec_epf_quad_wvar(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num = 1);


matrix objec_test(matrix par, dynprob &res1, data &dat1, matrix &randus, int nrun);

matrix simann(matrix x, bool max, double rt, double eps, int ns, int nt, int neps,
              int maxevl, matrix lb, matrix ub, matrix c, int iprint, int iseed1,
              double t, matrix vm, data &dat1, dynprob &res1, matrix mn, matrix W,
              matrix &momv, matrix &randus, string des, string direct,
              double (*objec) (matrix, matrix, matrix, dynprob&,
                      data&, matrix&, int, matrix&, int), double &fret, double &f_0, int thr_num = 1, int nopt = 0);

double objec_nm(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus,
                int nrun, matrix minb, matrix maxb, matrix &mn,
                    double (*objec) (matrix, matrix, matrix, dynprob&, data&, matrix&, int, matrix&, int), int thr_num = 1);



matrix out_auxiliary(matrix par, dynprob &res1, data &dat1, matrix &randus, int nrun, int thr_num,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int));

void out_auxiliaries(matrix par, dynprob &res1, data &dat1, matrix &randus, int nrun, int thr_num,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), matrix &report, matrix &report2);

matrix jacobian(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix mn,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval);

matrix variance_ineff(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix mn,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, matrix omega, int length);

matrix variance_eff(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix mn,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, int length);

matrix variance_mom(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix mn,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, matrix omega);

matrix variance_mom_outofsample(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix mn1, matrix mn2,
matrix (*auxiliary_est1)(dynprob &, int, int, data &, matrix&, int),
  matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, matrix omega1, matrix omega2, matrix omega12);

matrix jstat_ineff(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix m, matrix mn,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, matrix omega);

matrix jstat_eff(matrix m, matrix mn, data &dat1, matrix w);

matrix jstat_outofsample(matrix par, dynprob &res1, data &dat1, matrix &randus, matrix m1, matrix m2, matrix mn1, matrix mn2,
matrix (*auxiliary_est1)(dynprob &, int, int, data &, matrix&, int),
  matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix w, matrix omega1, matrix omega2, matrix omega12);


matrix test_stats(matrix par, matrix par_true, matrix v_par, matrix m, matrix m_true, matrix v_m, matrix &par_t, matrix &m_t);

void jacobians(matrix par, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0 = 1e-4, int l=-1);

void simple_grad(matrix par, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0 = 1e-4, int l = -1);

void jacobians_f(matrix par, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0 = 1e-4);

void jacobians_5pt(matrix par, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0 = 1e-4);

void full_inf(matrix par, matrix w, matrix omega1, matrix omega2, matrix omega12, int thr_num,
                dynprob res1, data dat1, matrix randus, int lengths, bool eff, matrix m, matrix m2,
                matrix (*aux1)(dynprob &, int, int, data &, matrix&, int),
                matrix (*aux2)(dynprob &, int, int, data &, matrix&, int),
                matrix &mn, matrix &mn2, matrix &var, matrix &var2, matrix &var3, matrix &jstat, matrix &jstat2);

matrix jstats(matrix m, matrix mn, matrix m2, matrix mn2, matrix vcvm, matrix v5, matrix v6, bool eff, matrix &jstat, matrix &jstats2, data dat1);

matrix test_stats(matrix par, matrix par_true, matrix v_par, matrix m, matrix m_true, matrix v_m, matrix &par_t, matrix &m_t);

#endif
