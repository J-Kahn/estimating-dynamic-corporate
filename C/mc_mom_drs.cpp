#include "matrix.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "epf.h"
#include "opt.h"
#include "regress.h"
#include <sstream>
#include <omp.h>
#include <sstream>
#include <mpi.h>
#include <limits>
#include "neldmead.h"

#define _USE_MATH_DEFINES
typedef std::numeric_limits< double > dbl;



int main(int argc, char **argv){


	// Start up the MPI part
  MPI_Init(&argc, &argv);

  // OpenMP threads
  int thr_num = 14;

	int id;
	int nth;

	// Start at some number in the Monte-Carlo trial data, do some number of trials
	int start = 0;
	int ntrials = 100;

	// Set initial point to the actual value used
  matrix init_points = readcsv("../Results/Estim/par_quad_drs_toni_rr3.csv",7,1);


	// Set up MPI grid
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &nth);

	// How many are we doing for each MPI?
	int nper = (ntrials - start) / nth;


	int nrun = 2;

	cout << "HERE TO START " << endl;

	// Read in the MC trial moments
	matrix moms = readcsv("../Results/MC/trial_data_rr_drs_mom.csv", 100, 8);
	matrix vmoms = readcsv("../Results/MC/trial_data_rr_drs_wmom.csv", 100, 8 * 8);

	cout << "HERE TO START 1" << endl;
	matrix lengths = readcsv("../Data/lengths.csv",2,1);

	int nfirms  = (int) lengths(1);
	int   nobs  = (int) lengths(0);
	int nobsper = nobs / nfirms;
  nfirms = 75000 / nobsper;

	cout << "HERE TO START 22" << endl;
	// Make the matrix of random draws
	matrix randus(nobsper * nfirms + 100 * nfirms + 3 * nfirms, 10);
	int seedn = 1158;
	srand(seedn);
	randus.randu();

  // Set up data and dynamic problem
  data dat1 = data(100, 100 + nobsper, nfirms, 60);
dynprob res1 = dynprob(15, 90, 15, 90, 9);
dat1.na2 = 60;

// Initialize interest rate
res1.pars[0]=1/(1+0.02 * (1 + 0.1));


// Initialize taxes on corporate profits
res1.pars[14] = 0.1;
res1.pars[21] = 0.0;

// Set up the order of estimation
res1.porder[0] = 2;
res1.porder[1] = 4;
res1.porder[2] = 9;
res1.porder[3] = 11;
res1.porder[4] = 1;
res1.porder[5] = 22;
res1.porder[6] = 23;

	// Bounds on parameters
	matrix minb(7), maxb(7);
	minb(0)=0.14,  maxb(0)=0.15;
	minb(1)=0.005, maxb(1)=0.095;
	minb(2)=0.45,  maxb(2)=0.55;
	minb(3)=0.15,  maxb(3)=0.25;
	minb(4)=0.5,  maxb(4)=0.6;
	minb(5)=0.25,  maxb(5)=0.35;
  minb(6)=0.45,  maxb(6)=0.55;
	matrix report;
	matrix fepf(8);


	// Set up parameters for simulated annealing
	matrix cb(7), par1, vm1(7);
	cb = cb + 2;
	vm1 = (maxb - minb)/2.0;
	double f, f_0;

	MPI_Barrier(MPI_COMM_WORLD);

	// Set up file names to save to
	ostringstream convert, append;
string res, filename, filername, resr;
	matrix vcvm, vcvmi;

	append << "_" << id;
	resr = append.str();


	filername = "../Results/MC/Trials/summary_mom_drs_rr10" + resr;

  ofstream filer;


	filer.open(filername.c_str(),std::ios_base::app);
	filer.precision(dbl::digits10);
  matrix var, var2, var3, jstat, jstat2, mn2, m, mn(moms.ncols());


            matrix dhmat(6);
            dhmat(0) = 1e-1;
            dhmat(1) = 1e-2;
            dhmat(2) = 1e-3;
            dhmat(3) = 1e-4;
            dhmat(4) = 1e-5;
            dhmat(5) = 1e-6;

	// Start running MC estimations
	for(int run = id * nper + start; run < (id +1) * nper + start; run++){
	// Run estimation


  		convert.str("");
  		convert.clear();
  		convert << id << "_" << run;
  		res = convert.str();
  		filename = "mom_drs_rr10" + res;

  		// Grab moments and W
  		m = moms.row(run);
  		vcvm = vmoms.row(run);
  		vcvm.reshape(moms.ncols(), moms.ncols());
      cout << "Moment " << id << ", " << run << + " ";
      for(int i = 0; i < moms.ncols(); i++){
        cout << m(i) << ", ";
      }
      cout << endl;

      cout << "Moment var "<< id << ", " << run << + " ";
      for(int i = 0; i < moms.ncols(); i++){
        for(int j = 0; j < moms.ncols(); j++){
          cout << vcvm(i, j) << ", ";
        }
      }

      cout << endl;
      res1.reset();
    matrix results = simann(init_points, false, 0.85, 1e-6, 20, 5, -4, 150000, minb, maxb, cb, 4, 100, 0.005,
	                       vm1, dat1, res1, m, vcvm, mn, randus, filename, "../Results/MC/Trials/", objec_mom, f, f_0, thr_num, 4);
      matrix jac, jac2, jac5, jac52;

      res1.wipe_long();
      matrix res_old = results;
    double fold = f;

      res1.wipe_long();

      cout << "DONE WITH RESULTS " << id << ", " << run << endl;
      results.print("../Results/Test/resul_" + filename + ".csv");
      parameter_update(results, res1);
      cout << "DONE WITH PARAMETER UPDATE " << id << ", " << run << endl;


      res1.vfi(1e-5,1e+16,3000,0,0,0,thr_num);

      cout << "DONE WITH VFI " << id << ", " << run << endl;


      res1.policy();

      cout << "DONE WITH POLICY " << id << ", " << run << endl;

      report = report_all_out_est(res1, dat1.nsim, dat1.nbr, dat1, randus, 10);

      cout << "DONE WITH REPORT " << id << ", " << run << endl;
      report.print("../Results/Test/repor_" + filename + ".csv");


      jacobians(results, res1, dat1, randus, mom_tm_est, epf_quad_est, 10, thr_num, 4, 4, jac, jac2, 1e-6, 4);
      cout << "DONE WITH JACOBIANS " << id << ", " << run << endl;
      jac.print("../Results/Test/jac1_" + filename + ".csv");
      jac2.print("../Results/Test/jac2_" + filename + ".csv");


	    // Save to a file

	    filer << run << ", " << f << ", " << f_0 << ", ";

    	for(int i = 0; i < results.nrows(); i++){
	      filer << results(i) << ", ";
    	}

	for(int i = 0; i < report.nrows(); i++){
	      filer << report(i) << ", ";
    	}

    	for(int i = 0; i < jac.nrows(); i++){
    	  for(int j = 0; j < jac.ncols(); j++){
    	    filer << jac(i,j) << ", ";
    	  }
    	}

    	for(int i = 0; i < jac2.nrows(); i++){
    	  for(int j = 0; j < jac2.ncols(); j++){
    	    filer << jac2(i,j) << ", ";
    	  }
    	}

    	for(int i = 0; i < mn.nrows(); i++){
	      filer << mn(i) << ", ";
    	}

      cout << "DONE WITH PRINTING " << id << ", " << run << endl;

    	filer << endl;


  }

	filer.close();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
