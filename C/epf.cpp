#include "matrix.h"
#include <cmath>
#include <omp.h>
#include "epf.h"
#include "regress.h"
#include "data.h"

#include <iostream>
#include <fstream>
#include <sstream>

//#include "UnivariateDensityDerivative.h"
#define _USE_MATH_DEFINES
using namespace std;




/*------------------------------------------------


            EPF/moment calculation


-------------------------------------------------*/

/*
    Polynomial EPF
                    */
matrix epf_quad(data &dat, int degree){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_quad(bt, rt, degree);
    matrix beta_i    = regress(it, x_regress);
    matrix beta_d    = regress(dt, x_regress);
    matrix beta_c    = regress(ct, x_regress);

    matrix beta(beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        beta(i)          = beta_i(i);
    }

    for(int i = 0; i < beta_d.nrows(); i++){
        beta(i + beta_i.nrows()) = beta_d(i);
    }

    for(int i = 0; i < beta_c.nrows(); i++){
        beta(i + beta_i.nrows() + beta_d.nrows()) = beta_c(i);
    }

    return beta;
}

/*
    Polynomial EPF, transformed state
                    */
matrix epf_quad_trans(data &dat, int degree){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("nwt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_quad(bt, rt, degree);
    matrix beta_i    = regress(it, x_regress);
    matrix beta_d    = regress(dt, x_regress);
    matrix beta_c    = regress(ct, x_regress);

    matrix beta(beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        beta(i)          = beta_i(i);
    }

    for(int i = 0; i < beta_d.nrows(); i++){
        beta(i + beta_i.nrows()) = beta_d(i);
    }

    for(int i = 0; i < beta_c.nrows(); i++){
        beta(i + beta_i.nrows() + beta_d.nrows()) = beta_c(i);
    }

    return beta;
}

/*
    Polynomial EPF, cubic form
                    */
matrix epf_cube(data &dat, int degree){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_cube(bt, rt, degree);
    matrix beta_i    = regress(it, x_regress);
    matrix beta_d    = regress(dt, x_regress);
    matrix beta_c    = regress(ct, x_regress);

    matrix beta(beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        beta(i)          = beta_i(i);
    }

    for(int i = 0; i < beta_d.nrows(); i++){
        beta(i + beta_i.nrows()) = beta_d(i);
    }

    for(int i = 0; i < beta_c.nrows(); i++){
        beta(i + beta_i.nrows() + beta_d.nrows()) = beta_c(i);
    }

    return beta;
}


/*
    Default moments
                    */
matrix mom_tm(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix moments(8);

    moments(2)  = mean(it);
    moments(3)  = var(it);
    moments(0)  = mean(bt);
    moments(1)  = var(bt);
    moments(4)  = mean(ct);
    moments(5)  = var(ct);
    moments(6)  = mean(rt);
    moments(7)  = var(rt);

    return moments;
}

/*
    Moments with covariance matrix
                    */
matrix mom_tm_cv(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix moments(11);

    moments(2)  = mean(it);
    moments(3)  = var(it);
    moments(0)  = mean(bt);
    moments(1)  = var(bt);
    moments(4)  = mean(ct);
    moments(5)  = var(ct);
    moments(6)  = mean(rt);
    moments(7)  = var(rt);
    moments(8)  = cov_vec(rt, bt);
    moments(9)  = cov_vec(rt, it);
    moments(10) = cov_vec(rt, dt);

    return moments;
}


/*------------------------------------------------


            Multiple simulations for estimation, takes
            average over simulations for SMM objective function.


-------------------------------------------------*/



matrix epf_quad_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){

    matrix ri = randus.col(0);
    sim(valstor,  ri, ret);
    matrix beta(18), betat(18);

    beta = epf_quad(ret, 2)/((double) ndat);

    for(int i = 1; i < ndat; i++){

        ri = randus.col(i);
        sim(valstor,  ri, ret);
        betat = epf_quad(ret, 2)/((double) ndat);
        if(beta.nrows() == betat.nrows()){
            beta = beta + betat;
        }
        else{
            beta(0) = INFINITY;
            cout << "Caught error for single" << endl;
        }

    }

    return beta;
}

matrix epf_quad_trans_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){

    matrix ri = randus.col(0);
    sim(valstor,  ri, ret);
    matrix beta(18), betat(18);

    beta = epf_quad_trans(ret, 2)/((double) ndat);

    for(int i = 1; i < ndat; i++){

        ri = randus.col(i);
        sim(valstor,  ri, ret);
        betat = epf_quad_trans(ret, 2)/((double) ndat);
        if(beta.nrows() == betat.nrows()){
            beta = beta + betat;
        }
        else{
            beta(0) = INFINITY;
            cout << "Caught error for single" << endl;
        }

    }

    return beta;
}

matrix epf_cube_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){

    matrix ri = randus.col(0);
    sim(valstor,  ri, ret);
    matrix beta(30), betat(30);

    beta = epf_cube(ret, 2)/((double) ndat);

    for(int i = 1; i < ndat; i++){

        ri = randus.col(i);
        sim(valstor,  ri, ret);
        betat = epf_cube(ret, 2)/((double) ndat);
        if(beta.nrows() == betat.nrows()){
            beta = beta + betat;
        }
        else{
            beta(0) = INFINITY;
            cout << "Caught error for single" << endl;
        }

    }

    return beta;
}


matrix report_all_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){

    matrix beta, betat;

    matrix ri = randus.col(0);
    sim(valstor,  ri, ret);
    beta = report_all_mc(valstor, ret.nsim, ret.nbr, ret, ri) / ((double) ndat);

    for(int i = 1; i < ndat; i++){

        ri = randus.col(i);
        sim(valstor,  ri, ret);
        betat = report_all_mc(valstor, ret.nsim, ret.nbr, ret, ri) / ((double) ndat);

        if(beta.nrows() == betat.nrows()){
            beta = beta + betat;
        }
        else{
            beta(0) = INFINITY;
            cout << "Caught error for single" << endl;
        }

    }

    return beta;
}

matrix report_all_out_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){

    matrix beta, betat;

    matrix ri = randus.col(0);
    sim(valstor,  ri, ret);
    beta = report_all_mc_out(valstor, ret.nsim, ret.nbr, ret, ri) / ((double) ndat);

    for(int i = 1; i < ndat; i++){

        ri = randus.col(i);
        sim(valstor,  ri, ret);
        betat = report_all_mc_out(valstor, ret.nsim, ret.nbr, ret, ri) / ((double) ndat);

        if(beta.nrows() == betat.nrows()){
            beta = beta + betat;
        }
        else{
            beta(0) = INFINITY;
            cout << "Caught error for single" << endl;
        }

    }

    return beta;
}



matrix mom_tm_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){
    matrix beta(8);
    matrix ri = randus.col(0);
    for(int i = 0; i < ndat; i++){
        ri = randus.col(i);
        sim(valstor,  ri, ret);
        beta = beta + mom_tm(ret)/((double) ndat);
    }
    return beta;
}



matrix mom_tm_cv_est(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus, int ndat){
    matrix beta(11);
    matrix ri = randus.col(0);
    for(int i = 0; i < ndat; i++){
        ri = randus.col(i);
        sim(valstor,  ri, ret);
        beta = beta + mom_tm_cv(ret)/((double) ndat);
    }
    return beta;
}


/*------------------------------------------------


            Influence functions for each estimator,
            for inputs into variance calculation.


-------------------------------------------------*/

matrix mom_tmINFL(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

	matrix moments(it.nrows(),8);

    matrix im = it - mean(it);
    matrix bm = bt - mean(bt);
    matrix cm = ct - mean(ct);
    matrix rm = rt - mean(rt);


    moments.col_sub(im,2);
    moments.col_sub((im & im) - var(it),3);
    moments.col_sub(bm,0);
    moments.col_sub((bm & bm) - var(bt),1);
    moments.col_sub(cm,4);
    moments.col_sub((cm & cm) - var(ct),5);
    moments.col_sub(rm, 6);
    moments.col_sub((rm & rm) - var(rt), 7);

	return moments;
}

matrix mom_tm_cvINFL(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix moments(it.nrows(),11);

    matrix im = it - mean(it);
    matrix bm = bt - mean(bt);
    matrix cm = ct - mean(ct);
    matrix dm = dt - mean(dt);
    matrix rm = rt - mean(rt);


    moments.col_sub(im,2);
    moments.col_sub((im & im) - var(it),3);
    moments.col_sub(bm,0);
    moments.col_sub((bm & bm) - var(bt),1);
    moments.col_sub(cm,4);
    moments.col_sub((cm & cm) - var(ct),5);
    moments.col_sub(rm, 6);
    moments.col_sub((rm & rm) - var(rt), 7);
    moments.col_sub((rm & bm) - cov_vec(rt, bt), 8);
    moments.col_sub((rm & im) - cov_vec(rt, it), 9);
    moments.col_sub((rm & dm) - cov_vec(rt,dt), 10);

    return moments;
}

matrix mom_tmINFL_nodemean(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix btdm = dat.get_var("bt"),
    rtdm        = dat.get_var("rt"),
    itdm        = dat.get_var("it"),
    dtdm        = dat.get_var("dt"),
    ctdm        = dat.get_var("ct");

    matrix moments(it.nrows(), 8);

    matrix im = it - mean(it);
    matrix bm = bt - mean(bt);
    matrix cm = ct - mean(ct);
    matrix rm = rt - mean(rt);

    matrix imdm = itdm - mean(itdm);
    matrix bmdm = btdm - mean(btdm);
    matrix cmdm = ctdm - mean(ctdm);
    matrix rmdm = rtdm - mean(rtdm);

    moments.col_sub(imdm,2);
    moments.col_sub((im & im) - var(it),3);
    moments.col_sub(bmdm,0);
    moments.col_sub((bm & bm) - var(bt),1);
    moments.col_sub(cmdm,4);
    moments.col_sub((cm & cm) - var(ct),5);
    moments.col_sub(rmdm, 6);
    moments.col_sub((rm & rm) - var(rt), 7);

    return moments;
}

/*
Moments with covariances
*/
matrix mom_tm_cvINFL_nodemean(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix btdm = dat.get_var("bt"),
    rtdm        = dat.get_var("rt"),
    itdm        = dat.get_var("it"),
    dtdm        = dat.get_var("dt"),
    ctdm        = dat.get_var("ct");

    matrix moments(it.nrows(), 11);

    matrix im = it - mean(it);
    matrix bm = bt - mean(bt);
    matrix dm = dt - mean(dt);
    matrix cm = ct - mean(ct);
    matrix rm = rt - mean(rt);

    matrix imdm = itdm - mean(itdm);
    matrix bmdm = btdm - mean(btdm);
    matrix cmdm = ctdm - mean(ctdm);
    matrix rmdm = rtdm - mean(rtdm);

    moments.col_sub(imdm,2);
    moments.col_sub((im & im) - var(it),3);
    moments.col_sub(bmdm,0);
    moments.col_sub((bm & bm) - var(bt),1);
    moments.col_sub(cmdm,4);
    moments.col_sub((cm & cm) - var(ct),5);
    moments.col_sub(rmdm, 6);
    moments.col_sub((rm & rm) - var(rt), 7);
    moments.col_sub((rm & bm) - cov_vec(rt, bt), 8);
    moments.col_sub((rm & im) - cov_vec(rt, it), 9);
    moments.col_sub((rm & dm) - cov_vec(rt, dt), 10);

    return moments;
}

/*
Quadratic EPF
*/
matrix epf_quadINFL(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_quad(bt, rt, 2);

    matrix i = it;
    matrix d = dt;
    matrix c = ct;


    matrix beta_i = solve_sym(x_regress.t()*x_regress,x_regress.t()*i);
    matrix beta_d = solve_sym(x_regress.t()*x_regress,x_regress.t()*d);
    matrix beta_c = solve_sym(x_regress.t()*x_regress,x_regress.t()*c);
    matrix res_i = i - x_regress * beta_i;
    matrix res_d = d - x_regress * beta_d;
    matrix res_c = c - x_regress * beta_c;

    matrix xr = solve_sym(x_regress.t() * x_regress, x_regress.t()).t()*((double) it.nrows());

    matrix infl(it.nrows(),beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_i,i);
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_d,i+beta_i.nrows());
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_c,i+2*beta_c.nrows());
    }


	return infl;
}

/*
Quadratic EPF with transformed state variables
*/
matrix epf_quad_transINFL(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("nwt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_quad(bt, rt, 2);

    matrix i = it;
    matrix d = dt;
    matrix c = ct;


    matrix beta_i = solve_sym(x_regress.t()*x_regress,x_regress.t()*i);
    matrix beta_d = solve_sym(x_regress.t()*x_regress,x_regress.t()*d);
    matrix beta_c = solve_sym(x_regress.t()*x_regress,x_regress.t()*c);
    matrix res_i = i - x_regress * beta_i;
    matrix res_d = d - x_regress * beta_d;
    matrix res_c = c - x_regress * beta_c;

    matrix xr = solve_sym(x_regress.t() * x_regress, x_regress.t()).t()*((double) it.nrows());

    matrix infl(it.nrows(),beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_i,i);
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_d,i+beta_i.nrows());
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_c,i+2*beta_c.nrows());
    }


	return infl;
}


/*
EPF with cubic polynomial.
*/
matrix epf_cubeINFL(data dat){

    matrix bt = dat.get_var("bt"),
    rt        = dat.get_var("rt"),
    it        = dat.get_var("it"),
    dt        = dat.get_var("dt"),
    ct        = dat.get_var("ct");

    matrix x_regress = tensor_cube(bt, rt, 2);

    matrix i = it;
    matrix d = dt;
    matrix c = ct;


    matrix beta_i = solve_sym(x_regress.t()*x_regress,x_regress.t()*i);
    matrix beta_d = solve_sym(x_regress.t()*x_regress,x_regress.t()*d);
    matrix beta_c = solve_sym(x_regress.t()*x_regress,x_regress.t()*c);
    matrix res_i = i - x_regress * beta_i;
    matrix res_d = d - x_regress * beta_d;
    matrix res_c = c - x_regress * beta_c;

    matrix xr = solve_sym(x_regress.t() * x_regress, x_regress.t()).t()*((double) it.nrows());

    matrix infl(it.nrows(),beta_i.nrows() + beta_d.nrows() + beta_c.nrows());

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_i,i);
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_d,i+beta_i.nrows());
    }

    for(int i = 0; i < beta_i.nrows(); i++){
        infl.col_sub(xr.col(i) & res_c,i+2*beta_c.nrows());
    }


	return infl;
}

/*------------------------------------------------


            Tools for turning influence functions
            into covariance matrices for moments.


-------------------------------------------------*/


/*
This computes variance matrices without clustering. It's not used anymore.
*/
matrix var_infl(data dat, matrix (*infl) (data)){
    matrix infli = infl(dat);
    matrix v = infli.cross();

    matrix w = v / (double) infli.nrows();

    return w;
}


/*
This clusters influence functions.
*/
matrix clust(matrix &infls, int nf, int nsim, int nbr){
    matrix inflc(nf, infls.ncols());
    #pragma omp parallel for
    for(int i = 0; i < nf; i++){
        for(int j = 0; j < infls.ncols(); j++){
            for(int t = 0; t < nsim
                - nbr; t++){
                inflc(i, j) += infls(i * (nsim-nbr) + t, j);
            }
        }
    }
    return inflc;
}

/*
This computes clustered variance matrices for one set of parameters.
*/
matrix var_infl_clust(data dat, matrix (*infl) (data)){
    matrix infli = infl(dat);
    matrix inflc  = clust(infli, dat.nf, dat.nsim, dat.nbr);
    matrix v = inflc.cross();
    matrix w = v / (double) infli.nrows();

    return w;
}

/*
This computes clustered covariance matrices for two sets of parameters
*/
matrix cov_infl_clust(data dat, matrix (*infl) (data), matrix (*infl2) (data)){
    matrix infli = infl(dat);
    matrix infli2 = infl2(dat);
    matrix inflc2 = clust(infli2, dat.nf, dat.nsim, dat.nbr);
    matrix inflc  = clust(infli, dat.nf, dat.nsim, dat.nbr);
    matrix v = inflc % inflc2;
    matrix w = v / (double) infli.nrows();

    return w;
}

/*------------------------------------------------


            Tools for producing a bunch of auxiliary
            estimators for a single trial, used in
            simulation.


-------------------------------------------------*/


/*
This prints off the moments and their covariance matrices for the data for Monte Carlo trials.
*/
matrix report_all_mc(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus){
    sim(valstor,  randus, ret);
    matrix beta_quad, beta_quad_wvar, beta_quad_wvarex, v_quad_wvarex2, beta_mom_kt, beta_mom_hp, beta_mom, betat, v_mom, v_mom_kt, v_mom_hp, v_quad, v_quad_wvar, v_quad_wvarex, beta_cube, beta_mom_cv, v_mom_cv, v_cube, v_momcube, v_mom_cvquad, v_mom_cv2, v_cube2, v_quad_trans, v_quad_trans2, beta_quad_trans, v_momquad_trans;
    matrix v_quad2, v_quad_wvar2, v_mom2, v_mom_kt2, v_mom_hp2, v_momquad, v_mom_ktquad, v_mom_hpquad, v_momquad_wvar, v_mom_ktquad_wvar, v_mom_hpquad_wvar, v_momquad_wvarex, v_mom_ktquad_wvarex, v_mom_hpquad_wvarex;

    /*
      Calculate relevant estimators.
    */

    beta_quad = epf_quad(ret, 2);

    beta_cube = epf_cube(ret, 2);

    beta_quad_trans = epf_quad_trans(ret, 2);

    beta_mom_cv = mom_tm_cv(ret);

    v_mom = var_infl_clust(ret, mom_tmINFL_nodemean);

    v_mom_cv = var_infl_clust(ret, mom_tm_cvINFL_nodemean);


    v_quad = var_infl_clust(ret, epf_quadINFL);

    v_cube = var_infl_clust(ret, epf_cubeINFL);

    v_mom2 = var_infl(ret, mom_tmINFL);

    v_cube2 = var_infl(ret, epf_cubeINFL);

    v_mom_cv2 = var_infl(ret, mom_tm_cvINFL);

    v_quad2 = var_infl(ret, epf_quadINFL);

    v_quad_trans = var_infl_clust(ret, epf_quad_transINFL);

    v_quad_trans2 = var_infl(ret, epf_quad_transINFL);

    v_momquad = cov_infl_clust(ret, epf_quadINFL, mom_tmINFL);

    v_mom_cvquad = cov_infl_clust(ret, epf_quadINFL, mom_tm_cvINFL);

    v_momcube = cov_infl_clust(ret, epf_cubeINFL, mom_tmINFL);

    v_momquad_trans = cov_infl_clust(ret, epf_quad_transINFL, mom_tmINFL);

    int ntot = beta_quad_trans.nrows() +beta_mom_cv.nrows() + beta_cube.nrows() +beta_quad.nrows() + beta_mom.nrows() + beta_quad_wvar.nrows() + beta_quad_wvarex.nrows()
                  + v_quad_trans.nrows() * v_quad_trans.nrows() + v_mom.nrows() * v_mom.nrows() + v_quad.nrows() * v_quad.nrows() + v_cube.nrows() * v_cube.nrows() + v_mom_cv.nrows() * v_mom_cv.nrows() + v_quad_wvar.nrows() * v_quad_wvar.nrows()+ v_quad_wvar.nrows() * v_quad_wvar.nrows() + v_mom_kt.nrows() * v_mom_kt.nrows() + v_mom_hp.nrows() * v_mom_hp.nrows()
                 + v_quad_trans.nrows() * v_quad_trans.nrows() + v_mom.nrows() * v_mom.nrows() + v_quad.nrows() * v_quad.nrows() + v_quad_wvar.nrows() * v_quad_wvar.nrows()+ v_quad_wvar.nrows() * v_quad_wvar.nrows() + v_mom_kt.nrows() * v_mom_kt.nrows() + v_mom_hp.nrows() * v_mom_hp.nrows()
                  + v_quad_trans.nrows() * v_mom.nrows() + v_quad.nrows() * v_mom_hp.nrows()+ v_quad.nrows() * v_mom_kt.nrows()+ v_quad.nrows() * v_mom.nrows()+ beta_mom_kt.nrows()
                 + beta_mom_hp.nrows() + v_quad_wvar.nrows() * v_mom_hp.nrows()+ v_quad_wvar.nrows() * v_mom_kt.nrows()+ v_quad_wvar.nrows() * v_mom.nrows()+v_quad_wvar.nrows() * v_mom_hp.nrows()+ v_quad_wvar.nrows() * v_mom_kt.nrows()+ v_quad_wvar.nrows() * v_mom.nrows() + v_cube.nrows() * v_cube.nrows() + v_mom_cv.nrows() * v_mom_cv.nrows() + + v_cube.nrows() * v_mom.nrows() + v_mom_cv.nrows() * v_quad.nrows();

    matrix beta(ntot);

    /*
       Return estimators as one long vector
    */

    int start = 0;

    for(int i = 0; i < beta_mom.nrows(); i++){
        beta(i) = beta_mom(i);
    }

    start += beta_mom.nrows();


    for(int i = 0; i < beta_quad.nrows(); i++){
        beta(i + start) = beta_quad(i);
    }

    start += beta_quad.nrows();


    for(int i = 0; i < beta_cube.nrows(); i++){
        beta(i + start) = beta_cube(i);
    }

    start += beta_cube.nrows();

    for(int i = 0; i < beta_mom_cv.nrows(); i++){
        beta(i + start) = beta_mom_cv(i);
    }

    start += beta_mom_cv.nrows();

    for(int i = 0; i < beta_quad_trans.nrows(); i++){
        beta(i + start) = beta_quad_trans(i);
    }

    start += beta_quad_trans.nrows();

    // Clustered variances

    for(int i = 0; i < v_mom.nrows(); i++){
        for(int j = 0; j < v_mom.nrows(); j++){
            beta(i*v_mom.nrows() + j + start) = v_mom(i,j);
        }
    }

    start += v_mom.nrows() * v_mom.nrows();

    for(int i = 0; i < v_quad.nrows(); i++){
        for(int j = 0; j < v_quad.nrows(); j++){
            beta(i*v_quad.nrows() + j + start) = v_quad(i,j);
        }
    }

    start += v_quad.nrows() * v_quad.nrows();


    for(int i = 0; i < v_cube.nrows(); i++){
        for(int j = 0; j < v_cube.nrows(); j++){
            beta(i*v_cube.nrows() + j + start) = v_cube(i,j);
        }
    }

    start += v_cube.nrows() * v_cube.nrows();

    for(int i = 0; i < v_mom_cv.nrows(); i++){
        for(int j = 0; j < v_mom_cv.nrows(); j++){
            beta(i*v_mom_cv.nrows() + j + start) = v_mom_cv(i,j);
        }
    }

    start += v_mom_cv.nrows() * v_mom_cv.nrows();

    for(int i = 0; i < v_quad_trans.nrows(); i++){
        for(int j = 0; j < v_quad_trans.nrows(); j++){
            beta(i*v_quad_trans.nrows() + j + start) = v_quad_trans(i,j);
        }
    }

    start += v_quad_trans.nrows() * v_quad_trans.nrows();


    for(int i = 0; i < v_momquad.nrows(); i++){
        for(int j = 0; j < v_momquad.ncols(); j++){
            beta(i*v_momquad.ncols() + j + start) = v_momquad(i,j);
        }
    }

    start += v_momquad.nrows() * v_momquad.ncols();

    for(int i = 0; i < v_momcube.nrows(); i++){
        for(int j = 0; j < v_momcube.ncols(); j++){
            beta(i*v_momcube.ncols() + j + start) = v_momcube(i,j);
        }
    }

    start += v_momcube.nrows() * v_momcube.ncols();

    for(int i = 0; i < v_mom_cvquad.nrows(); i++){
        for(int j = 0; j < v_mom_cvquad.ncols(); j++){
            beta(i*v_mom_cvquad.ncols() + j + start) = v_mom_cvquad(i,j);
        }
    }

    start += v_mom_cvquad.nrows() * v_mom_cvquad.ncols();

    for(int i = 0; i < v_momquad_trans.nrows(); i++){
        for(int j = 0; j < v_momquad_trans.ncols(); j++){
            beta(i*v_momquad_trans.ncols() + j + start) = v_momquad_trans(i,j);
        }
    }

    start += v_momquad_trans.nrows() * v_momquad_trans.ncols();

    cout << ntot << endl;
    cout << start << endl;

    return beta;
}

/*
This prints out all the estimators for the results of the Monte Carlos.
*/
matrix report_all_mc_out(dynprob &valstor, int nsim, int nbr, data &ret, matrix&randus){

    sim(valstor,  randus, ret);

    matrix beta_quad, beta_quad_wvar, beta_wvar, beta_quad_wvarex, beta_mom_kt, beta_mom_hp, beta_mom, betat;

    beta_mom = mom_tm(ret);

    // These are just here because of deprecated estimators, they keep consistency with python code.
    beta_mom_kt = zeros(18,1);
    beta_mom_hp = zeros(18,1);

    beta_quad_wvar = zeros(21, 1);
    beta_quad_wvarex = zeros(21,1);


    beta_quad = epf_quad(ret, 2);

    int ntot = beta_quad.nrows() + beta_mom.nrows() + beta_quad_wvar.nrows() + beta_quad_wvarex.nrows() + beta_mom_kt.nrows()
                 + beta_mom_hp.nrows();

    matrix beta(ntot);

    int start = 0;

    for(int i = 0; i < beta_mom.nrows(); i++){
        beta(i) = beta_mom(i);
    }

    start += beta_mom.nrows();

    for(int i = 0; i < beta_mom_kt.nrows(); i++){
        beta(i+start) = beta_mom_kt(i);
    }

    start += beta_mom_kt.nrows();

    for(int i = 0; i < beta_mom_hp.nrows(); i++){
        beta(i+start) = beta_mom_hp(i);
    }

    start += beta_mom_hp.nrows();


    for(int i = 0; i < beta_quad.nrows(); i++){
        beta(i + start) = beta_quad(i);
    }

    start += beta_quad.nrows();

    for(int i = 0; i < beta_quad_wvar.nrows(); i++){
        beta(i + start) = beta_quad_wvar(i);
    }

    start += beta_quad_wvar.nrows();


    for(int i = 0; i < beta_quad_wvarex.nrows(); i++){
        beta(i + start) = beta_quad_wvarex(i);
    }

    return beta;
}
