#include "matrix.h"
#include <math.h>
#include "epf.h"
#include "solve.h"
#include "opt.h"
#include <omp.h>
#define _USE_MATH_DEFINES
using namespace std;

/*

This command is to update the parameters of the model every time we consider
a new potential parameter estimate.


*/

void parameter_update(matrix &par, dynprob &prob){

    /* For reference: beta             = pars(0),  // discount rate
					alpha            = pars(1),  // curvature
					delta            = pars(2),  // depreciation
					lambda0          = pars(3),  // fixed issuance cost
					lambda1          = pars(4),  // lambda
					lambda2          = pars(5),  // quadratic issuance cost
					phi              = pars(6),  // debt issuance cost
					flowcost         = pars(7),  // flow fixed cost
					r_p              = pars(8),  // premium to debt
					firesale         = pars(9),  // firesale value of capital
					fixc             = pars(10), // fixed adjustment cost
					smoothc          = pars(11), // smooth adjustment cost
					tau_c            = pars(12), // corporate tax rate
					tau_d            = pars(13); // distribution tax rate
					debtc_down       = pars(14); // downward debt issuance cost
					debtc_up         = pars(15); // upward debt issuance cost
					cashc            = pars(16); // cost for cash balances
					debtc_fixed_down = pars(17); // fixed debt adjustment cost
					debtc_fixed_up   = pars(18); // fixed debt adjustment cost
					debtc_adj_up     = pars(19); // downward debt adjustment cost
					debtc_adj_down   = pars(20); // upward debt adjustment cost
				            mu               = pars[21];
            rho              = pars[
            sigma            = pars[23];
					*/
    for(int i = 0; i < par.nrows(); i++){
    	prob.pars[prob.porder[i]] = prob.scale[prob.porder[i]] * par(i);
    }

    prob.lom_tauchen();

    prob.grid_setup();

    prob.interpolant();

    prob.flow();

}

void parameter_update(dynprob &prob){

    prob.lom_tauchen();

    prob.grid_setup();
    prob.interpolant();

    prob.flow();
}


/*

Objection for moments

Inputs:

par    : parameter vector.
m      : data moments/EPF.
W      : data variance covariance matrix.
res1   : dynprob object to store reslts of VFI.
dat1   : data object for simulation.
randus : random uniform matrix for simulations.
nrun   : number of simulations.
mn     : vector to store moments from model.

*/

double objec_mom(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num){
// Set pars: beta, alpha, delta, gamma, Gamma, r_f, r_p, lambda, ef
	try{

		parameter_update(par, res1);
	    res1.vfi(1e-5,1e+16,3000,-1,0,0,thr_num);
    	res1.policy();
            omp_set_num_threads(thr_num);
	    matrix did  = mom_tm_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
	    if(did.nrows() == mn.nrows()){
	        mn = did;
	    }
	    did = did - m;
	    matrix dist = (did.t())*W*did;

	    if(isnan(dist(0,0))){
		res1.reset();
		    dist(0,0)=INFINITY;
	        cout << "NAN encountered" << endl;
	    }
	    if(res1.iter>3000){
		res1.reset();
	        dist(0,0) = INFINITY;
	        cout << "VFI not solved" << endl;
	    }
	    return dist(0,0);
	}
	catch(int e){
	    res1.reset();
	    return INFINITY;
	}
}


double objec_mom_cv(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num){
// Set pars: beta, alpha, delta, gamma, Gamma, r_f, r_p, lambda, ef
	try{

		parameter_update(par, res1);
	    res1.vfi(1e-5,1e+16,3000,-1,0,0,thr_num);
    	res1.policy();

	    matrix did  = mom_tm_cv_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
	    if(did.nrows() == mn.nrows()){
	        mn = did;
	    }
	    did = did - m;
	    matrix dist = (did.t())*W*did;

	    if(isnan(dist(0,0))){
		res1.reset();
		    dist(0,0)=INFINITY;
	        cout << "NAN encountered" << endl;
	    }
	    if(res1.iter>3000){
		res1.reset();
	        dist(0,0) = INFINITY;
	        cout << "VFI not solved" << endl;
	    }
	    return dist(0,0);
	}
	catch(int e){
	    res1.reset();
	    return INFINITY;
	}
}


/*

Objection for polynomial EPF

Inputs:

par    : parameter vector.
m      : data moments/EPF.
W      : data variance covariance matrix.
res1   : dynprob object to store reslts of VFI.
dat1   : data object for simulation.
randus : random uniform matrix for simulations.
nrun   : number of simulations.
mn     : vector to store moments from model.

*/

double objec_epf_quad(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num){
    try{
		parameter_update(par, res1);
	    res1.vfi(1e-5,1e+16,3000,-1,0,0,thr_num);
	    res1.policy();


            omp_set_num_threads(thr_num);
	    matrix did  = epf_quad_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
	    matrix dist;
	    //cout << did.ncols() << " " << did.nrows() << endl;
	    if(did.nrows() == m.nrows()){
	        mn = did;
	        did = did -m;
	        dist = (did.t())*W*did;

	    }
	    else{
	    res1.wipe();
	        dist=matrix(1,1);
	        dist(0,0) = INFINITY;
	        cout << "IT HAPPENED AGAIN" << endl;
	    }
	    if(isnan(dist(0,0))){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "NAN encountered" << endl;
	    }
	    if(res1.iter>3000){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "VFI not solved" << endl;
	    }
	    return dist(0,0);
	}
	catch(int e){
	    res1.wipe();
	    return INFINITY;
	}
}

double objec_epf_cube(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num){
    try{

		parameter_update(par, res1);
	    res1.vfi(1e-5,1e+16,3000,-1,0,0,thr_num);
	    res1.policy();
omp_set_num_threads(thr_num);


	    matrix did  = epf_cube_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
	    matrix dist;
	    //cout << did.ncols() << " " << did.nrows() << endl;
	    if(did.nrows() == m.nrows()){
	        mn = did;
	        did = did -m;
	        dist = (did.t())*W*did;

	    }
	    else{
	    res1.wipe();
	        dist=matrix(1,1);
	        dist(0,0) = INFINITY;
	        cout << "IT HAPPENED AGAIN" << endl;
	    }
	    if(isnan(dist(0,0))){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "NAN encountered" << endl;
	    }
	    if(res1.iter>3000){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "VFI not solved" << endl;
	    }
	    return dist(0,0);
	}
	catch(int e){
	    res1.wipe();
	    return INFINITY;
	}
}


/*

Objection for polynomial EPF

Inputs:

par    : parameter vector.
m      : data moments/EPF.
W      : data variance covariance matrix.
res1   : dynprob object to store reslts of VFI.
dat1   : data object for simulation.
randus : random uniform matrix for simulations.
nrun   : number of simulations.
mn     : vector to store moments from model.

*/

double objec_epf_quad_trans(matrix par, matrix m, matrix W, dynprob &res1, data &dat1, matrix &randus, int nrun, matrix &mn, int thr_num){
    try{

		parameter_update(par, res1);
	    res1.vfi(1e-5,1e+16,3000,-1,0,0,thr_num);
	    res1.policy();

omp_set_num_threads(thr_num);

	    matrix did  = epf_quad_trans_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
	    matrix dist;
	    //cout << did.ncols() << " " << did.nrows() << endl;
	    if(did.nrows() == m.nrows()){
	        mn = did;
	        did = did -m;
	        dist = (did.t())*W*did;

	    }
	    else{
	    res1.wipe();
	        dist=matrix(1,1);
	        dist(0,0) = INFINITY;
	        cout << "IT HAPPENED AGAIN" << endl;
	    }
	    if(isnan(dist(0,0))){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "NAN encountered" << endl;
	    }
	    if(res1.iter>3000){
	    res1.wipe();
	        dist(0,0) = INFINITY;
	        cout << "VFI not solved" << endl;
	    }
	    return dist(0,0);
	}
	catch(int e){
	    res1.wipe();
	    return INFINITY;
	}
}


matrix out_auxiliary(matrix par, dynprob &res1, data &dat1, matrix &randus, int nrun, int thr_num,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int)){

        matrix randi = randus.col(0);

        parameter_update(par, res1);

        res1.vfi(1e-5,1e+16,3000,-1,0,0, thr_num);

        res1.policy();

	    	matrix report = auxiliary_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);

        return report;
}

void out_auxiliaries(matrix par, dynprob &res1, data &dat1, matrix &randus, int nrun, int thr_num,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), matrix &report, matrix &report2){

        matrix randi = randus.col(0);

        parameter_update(par, res1);

        res1.vfi(1e-5,1e+16,3000,-1,0,0, thr_num);

        res1.policy();
omp_set_num_threads(thr_num);

            dat1.wipe();
		    report  = auxiliary_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
            dat1.wipe();
		    report2 = auxiliary_est2(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
            dat1.wipe();

        /*if(report.nrows() < 8){
          report = matrix(8);
          for(int i = 0; i < 8; i++){
            report(i) = INFINITY;
          }
        }
        if(report2.nrows() < 18){
          report2 = matrix(18);
          for(int i = 0; i < 18; i++){
            report2(i) = INFINITY;
          }
        }*/
    res1.wipe_long();

}
