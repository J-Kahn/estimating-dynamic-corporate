#include "matrix.h"
#include <cmath>
#include "epf.h"
#include "opt.h"
#include "regress.h"
#include <limits>

#define _USE_MATH_DEFINES

using namespace std;
typedef std::numeric_limits< double > dbl;

double ranmar(){
    double r = ((double) rand())/((double) RAND_MAX);
    return r;
}

/*======================================================

A BUNCH OF PRINT FUNCTIONS FOR SIMULATED ANNEALING

========================================================*/

void prt1(){

    cout << "  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS " << endl <<
            "  (lb AND ub). execution terminated without any"  << endl <<
            "  optimization. respecify x, ub OR lb so that  "  << endl <<
            "  lb(i) < x(i) < ub(i), i = 1, n. " << endl;
}


void prt2(bool max, int n, matrix x, double f){

    cout << "Initial X" << endl;
    x.print();
    cout << endl;
    if(max){
        cout <<"  INITIAL F: " <<  f << endl;
    }
    else{

        cout <<"  INITIAL F: " <<  -f << endl;
    }
}


void prt3(bool max, int n, matrix xp, matrix x, double f){

    cout << "Current X" << endl;
    x.print();
    cout << endl;

    if(max){
      cout << "  CURRENT F: " << f << endl;
    }
    else{
      cout << "  CURRENT F: " << -f << endl;
    }
    cout << "Trial X" << endl;
    xp.print();
    cout << endl;
    cout << "  POINT REJECTED SINCE OUT OF BOUNDS" << endl;
}


void prt4(bool max, int n, matrix xp, matrix x, double fp, double f){
    cout << "Current X" << endl;
    x.print();
    cout << endl;
    if(max){
        cout  << "  CURRENT F: " << f << endl;
        cout << "Trial X" << endl;
        xp.print();
        cout << endl;
        cout  << "  RESULTING F: " << fp << endl;
    }
    else{
        cout  << "  CURRENT F: " << -f << endl;
        cout << "Trial X" << endl;
        xp.print();
        cout << endl;
        cout  << "  RESULTING F: " << -fp << endl;
    }
}


void prt5(){
    cout << "  TOO MANY FUNCTION EVALUATIONS; CONSIDER " << endl <<
               "  increasing maxevl OR eps, OR decreasing"  << endl <<
               "  nt OR rt. these results are likely TO be poor." << endl;
}


void prt6(bool max){
    if (max) {
        cout  << "  THOUGH LOWER, POINT ACCEPTED" << endl;
    }
    else{
        cout  << "  THOUGH HIGHER, POINT ACCEPTED"  << endl;
    }
}


void prt7(bool max){
    if (max) {
        cout  << "  LOWER POINT REJECTED" << endl;
    }
    else{
        cout  << "  HIGHER POINT REJECTED" << endl;
    }
}


void prt8(int n, matrix vm, matrix xopt, matrix x, string des, matrix report_all, string direct){
    cout  << " intermediate results after step length adjustment" << endl;
    cout << "NEW STEP LENGTH" << endl;
    vm.print();
    cout << endl;
    cout << "Current Optimal X" << endl;
    xopt.print();
    xopt.print(direct + "par_" + des + ".csv");
    report_all.print(direct + "report" + des + ".csv");
    cout << endl;
    cout << "Current X" << endl;
    x.print();
    cout << endl;
    vm.print(direct + "vm_" + des + ".csv");
}

void prt9(bool max, int n, double t, matrix xopt, matrix vm, double fopt, int nup, int ndown,
          int nrej, int lnobds, int nnew, int nobjecev, matrix momopt, int nmom, string des, string direct, double* vimmortal, int ns){
    
    int totmov = nup + ndown + nrej;
    vm.print();
    momopt.print(direct + "mom_" + des + ".csv");
    xopt.print(direct + "par_" + des + ".csv");
    doublePrint(vimmortal, "vim_" + des + ".csv", ns);
    cout << " intermediate results before next temperature reduction" << endl;
    cout << "  CURRENT TEMPERATURE:            " << t << endl;
    
    if(max){
        cout <<"  MAX FUNCTION VALUE SO FAR:  " << fopt << endl;
        cout << "  TOTAL MOVES:                " << totmov << endl;
        cout << "     UPHILL:                  " << nup << endl;
        cout << "     ACCEPTED DOWNHILL:       " << ndown << endl;
        cout << "     REJECTED DOWNHILL:       " << nrej << endl;
        cout << "  OUT OF BOUNDS TRIALS:       " << lnobds << endl;
        cout << "  NEW MAXIMA THIS TEMPERATURE:" << nnew << endl;
        cout << "  TOTAL NUMBER OF EVALUATIONS:" << nobjecev << endl;
    }
    else{
        cout <<"  MIN FUNCTION VALUE SO FAR:  " << -fopt << endl;
        cout <<"  TOTAL MOVES:                " << totmov << endl;
        cout <<"     DOWNHILL:                " <<  nup << endl;
        cout <<"     ACCEPTED UPHILL:         " <<  ndown << endl;
        cout <<"     REJECTED UPHILL:         " <<  nrej << endl;
        cout <<"  TRIALS OUT OF BOUNDS:       " <<  lnobds << endl;
        cout <<"  NEW MINIMA THIS TEMPERATURE:" <<  nnew << endl;
        cout <<"  TOTAL NUMBER OF EVALUATIONS:" << nobjecev << endl;
    }
    
    cout << "CURRENT OPTIMAL X" << endl;
    xopt.print();
    cout << endl;
    
    cout << "STEP LENGTH (VM)" << endl;
    vm.print();
}




void prt10(){
    cout << "  SA ACHIEVED TERMINATION CRITERIA. IER = 0. " << endl;
}

void prt110(double f, double fp, double fopt, matrix x, matrix xopt,matrix xd, matrix momv,
    matrix momopt, matrix momd, matrix m, double T, int k, int infer, matrix WL, string file){
      ofstream simfile;
       
      matrix md = WL.t() * (m - momopt);
      matrix md2 = WL.t() * (m - momv);
      matrix md3 = WL.t() * (m - momd);
      simfile.open(file.c_str());
      simfile.precision(dbl::digits10);
      simfile << "-------------------------------------- " << endl;
      simfile << "Energy: " << f << " comparison state: " << fp << " best: " << fopt << endl;
      simfile << "Trial: ";
      simfile << x(0);
      for(int i = 1; i < x.nrows(); i++){
         simfile << ", " << x(i);
      }
      simfile << endl;
      simfile << "Current: ";
      simfile << xd(0);
      for(int i = 1; i < xd.nrows(); i++){
         simfile << ", " << xd(i);
      }
      simfile << endl;
      simfile << "Best: ";
      simfile << xopt(0);
      for(int i = 1; i < xopt.nrows(); i++){
         simfile << ", " << xopt(i);
      }
      simfile << endl;
      simfile << endl;
      /*cout << "Simulated   Best        Actual ";
      for(int i = 0; i < momv.nrows(); i++){
        cout << setprecision(3) << momv(i) << "       " << momopt(i) << "       " << m(i) << endl;
      }
      cout << endl;*/
       
      simfile << "Trial: ";
      simfile << momv(0);
      for(int i = 1; i < momv.nrows(); i++){
         simfile << ", " << momv(i);
      }
      simfile << endl << endl;
      simfile << "Current: ";
      simfile << momd(0);
      for(int i = 1; i < momd.nrows(); i++){
         simfile << ", " << momd(i);
      }
      simfile << endl << endl;
      simfile << "Best: ";
      simfile << momopt(0);
      for(int i = 1; i < momopt.nrows(); i++){
         simfile << ", " << momopt(i);
      }
      simfile<< endl << endl;
      
      simfile << "Actual: ";
      simfile << m(0);
      for(int i = 1; i < m.nrows(); i++){
         simfile << ", " << m(i);
      }
      simfile << endl << endl;

      simfile << "W Trial: ";
      simfile << md2(0);
      for(int i = 1; i < md2.nrows(); i++){
         simfile << ", " << md2(i);
      }
      simfile << endl << endl;

      simfile << "W Current: ";
      simfile << md3(0);
      for(int i = 1; i < md3.nrows(); i++){
         simfile << ", " << md3(i);
      }
      simfile << endl << endl;
      simfile << "W Best: ";
      simfile << md(0);
      for(int i = 1; i < md.nrows(); i++){
         simfile << ", " << md(i);
      }
      simfile << endl << endl;

      simfile << "Temp: " << T << " evaluations " << k << " P " << exp((fp - f)/T) << " infer " << infer << endl;
      simfile << "-------------------------------------- " << endl << endl;
}


void prt11(double f, double fp, double fopt, matrix x, matrix xopt, matrix momv,
    matrix momopt, matrix m, double T, int k, int infer, string file, string file2, string des, string direct){
      cout << "-------------------------------------- " << endl;
      cout << "Energy: " << f << " comparison state: " << fp << " best: " << fopt << endl;
      cout << "Current: ";
      (x.t()).printv();
      cout << endl;
      cout << "Best: ";
      (xopt.t()).printv();
      cout << endl;
      xopt.print(direct + "par_" + des + ".csv");
      cout << endl;
      
      cout << "Simulated: ";
      (momv.t()).printv();
      cout << endl << endl;
      cout << "Best: ";
      (momopt.t()).printv();
      cout << endl << endl;
      cout << "Actual: ";
      (m.t()).printv();
      cout << endl << endl;

      cout << "Temp: " << T << " evaluations " << k << " P " << exp((fp - f)/T) << " infer " << infer << endl;
      cout << "-------------------------------------- " << endl << endl;

        ofstream simfile;
        simfile.open(file.c_str(), ios::app);
        simfile.precision(dbl::digits10);
        simfile << x(0);
        for(int i = 1; i < x.nrows(); i++){
            simfile << ", " << x(i);
        }
        simfile << ", " << fp << ", " << f << ", " << exp((fp - f)/T) << endl;
        simfile.close();


        ofstream simfile2;
        simfile2.open(file2.c_str(), ios::app);
        simfile2.precision(dbl::digits10);
        simfile2 << momopt(0);
        for(int i = 1; i < momopt.nrows(); i++){
            simfile2 << ", " << momopt(i);
        }
        simfile2 << endl;
        simfile2.close();
}

/*
Simulated annealing code. Finds minimum of objective function, objec, of the form in objec.cpp. Ported from Bill Goffe's Fortran code.

Arguments:

    N      : Number of variables in the function to be optimized. (INT)

    X      : The starting values for the variables of the function to be optimized. (DP(N))

    MAX    : Denotes whether the function should be maximized or minimized. A true value denotes maximization while a false value denotes minimization. (L)

    RT     : The temperature reduction factor.  The value suggested by Corana et al. is .85. See Goffe et al. for more advice. (DP)

    EPS    : Error tolerance for termination. (EP)

    NS     : Number of cycles.  After NS*N function evaluations, each element of VM is adjusted so that approximately half of all function evaluations are accepted.  The suggested value is 20. (INT)

    NT     : Number of iterations before temperature reduction.

    NEPS   : Number of final function values used to decide upon termination.  See EPS.  Suggested value is 4. (INT)

    MAXEVL : The maximum number of function evaluations.  If it is exceeded, IER = 1. (INT)

    LB     : The lower bound for the allowable solution variables. (DP(N))

    UB     : The upper bound for the allowable solution variables. (DP(N)) If the algorithm chooses X(I)

    C      : Vector that controls the step length adjustment.  The suggested value for all elements is 2.0. (DP(N))

    IPRINT : controls printing inside SA. (INT), 0 - 3, greater levels of printing.

    ISEED1 : The first seed for the random number generator RANMAR. 0 <= ISEED1 <= 31328. (INT)

    ISEED2 : The second seed for the random number generator RANMAR. 0 <= ISEED2 <= 30081.

    T      : On input, the initial temperature. See Goffe et al. for advice. On output, the final temperature. (DP)

    VM     : The step length vector. On input it should encompass the region of interest given the starting value X.

    dat1   : Data object to store simulations.

    res1   : Dynprob object to store solutions.

    mn     : Data moments vector.

    W      : Variance of data moments.

    momv   : Simulated moment vector.

    randus : Matrix of random numbers for simulations.

    des    : String name for saving files.

    object : Objective function to be minimized.

*/
matrix simann(matrix x, bool max, double rt, double eps, int ns, int nt, int neps,
              int maxevl, matrix lb, matrix ub, matrix c, int iprint, int iseed1,
              double t, matrix vm, data &dat1, dynprob &res1, matrix mn, matrix W,
              matrix &momv, matrix &randus, string des, string direct,
              double (*objec) (matrix, matrix, matrix, dynprob&,
                      data&, matrix&, int, matrix&, int), double &fret, double &f_0, int thr_num, int nopt
              ){

    matrix report_all(0);

    /*if(iprint > -1){
      report_all = matrix(48 + 18 + 8);
    }*/

    int nmom = momv.nrows();
    
    	matrix WL = chol(W);


    x.print();
    lb.print();
    ub.print();
    int nrun = randus.ncols();

    // We need this to access the file names
    int n = x.nrows();
    int infer = 0;

    if(nopt <= 0){
      nopt = n;
    }
   
    //  Type all external variables.
    matrix xopt(n);
    double fopt, diff;
    
    //  Type all internal variables.
    double f, fp=0, p, pp, ratio, fpl;
    matrix xp(n), fstar(neps);
    int nup, ndown, nrej, nnew, lnobds, h, i, j, m;
    imat nacp(n);
    bool quit;
    matrix vgopt(nmom), momopt(nmom), momd(nmom);

    srand(iseed1);

    //  Set initial values.
    int nacc = 0, nobds = 0, nobjecev = 0, ier = 99;
    
    for(int i = 0; i < n; i++){
        xopt(i) = x(i);
        nacp(i) = 0;
    }
    
    for(int i = 0; i < neps; i++){
        fstar(i) = 1e+15;
    }
    
    //  if the initial temperature is not positive, notify the user and
    //  return to the ing routine.
    if (t <= 0.0) {
        cout  << "  THE INITIAL TEMPERATURE IS NOT POSITIVE. " << endl <<
                 "  reset the variable t. " << endl;
        ier = 3;
        goto end;
    }

    //  if the initial value is out of bounds, notify the user and return
    //  to the ing routine.
    for(int i = 0; i < n; i++){
        if ((x(i) > ub(i)) || (x(i) < lb(i))) {
            prt1();
            ier = 2;
            goto end;
        }
    }

    //  Evaluate the function with input X and return value as F.
    fopt = 100000000;

    f = objec(x, mn, W, res1, dat1, randus, nrun, momv, thr_num); // OBJECTIVE!!!
	cout << "Objective" << endl;
    momopt = momv;
    momd = momv;
	cout << "Momv " << endl;
    //  if the function is to be minimized, switch the sign of the function.
    //  Note that all intermediate and final output switches the sign back
    //  to eliminate any possible confusion for the user.
    if(! max) f = -f;
    nobjecev++;
    fopt = f;
    f_0 = f;
    fstar(0) = f;
    mn.print(direct + "mom_init_" + des + ".csv");
    if(iprint >= 1)   prt2(max, n, x, f);

    //  Start the main loop. Note that it terminates if (i) the algorithm
    //  succ}}esfully optimizes the function or (ii) there are too many
    //  function evaluations (more than MAXEVL).
    max_loop:

    nup = 0, nrej = 0, nnew = 0, ndown = 0, lnobds = 0;
    
    for(int m = 0; m < nt; m++){
        for(int j = 0; j < ns; j++){
            for(int h = 0; h < nopt; h++){

                // PROBLEM AREA
                
                //  Generate XP, the trial value of X. Note use of VM to choose XP.
                for(int i = 0; i < n; i++){
                    if (i == h) {
                        xp(i) = x(i) + (ranmar()*2.0 - 1.0) * vm(i);
                    }
                    else{
                        xp(i) = x(i);
                    }

                    //  if XP is out of bounds, select a point in bounds for the trial.
                    if((xp(i) < lb(i)) || (xp(i) > ub(i))) {
                        xp(i) = lb(i) + (ub(i) - lb(i))*ranmar();
                        lnobds++;
                        nobds++;
                        if(iprint >= 3)   prt3(max, n, xp, x, f);
                    }
                }

                // END PROBLEM AREA

                fpl = fp;

                //  Evaluate the function with the trial point XP and return as FP.
                fp = objec(xp, mn, W, res1, dat1, randus, nrun, momv, thr_num);
                

		if(! max) fp = -fp;
                nobjecev++;
                if(iprint >= 3)   prt4(max, n, xp, x, fp, f);

                // SHOULD BE UNNECESSARY
                if(fp == INFINITY || fp == -INFINITY){
                    infer++;
                }
                else{
                    infer = 0;
                }
                if(infer > 10){
                    f = fopt;
                    x = xopt;
                    momd = momopt;
                    infer = 0;
                }
                // eND UNNECESSARY
                if(iprint >= 4){
                    prt11(f, fp,  fopt, xp, xopt, momv, momopt, mn, t, nobjecev, infer, direct + "simann_track_" + des + ".csv", direct + "simann_mom_" + des + ".csv", des, direct);
                }
                if(iprint >= -1){
                  if(nobjecev % 1000 == 0){
                    prt110(f, fp,  fopt, xp, xopt, x, momv, momopt, momd, mn, t, nobjecev, infer, WL, direct + "track_" + des + ".csv");
                  }
                }
                //  if too many function evaluations occur, terminate the algorithm.
                if(nobjecev >= maxevl) {
                    prt5();
                    if (! max) fopt = -fopt;
                    ier = 1;
                    goto end;
                }

                //  Accept the new point if the function value increases.
                if(fp >= f) {
                    if(iprint >= 3) {
                      cout  << "  POINT ACCEPTED" << endl;
                    }
                    x = xp;
	            momd = momv;
                    f = fp;
                    nacc++;
                    nacp(h) = nacp(h) + 1;
                    nup++;
                    
                    //  if greater than any other point, record as new optimum.
                    if (fp > fopt) {

                        if(report_all.nrows() > 0){
                          report_all = report_all_est(res1, dat1.nsim, dat1.nbr, dat1, randus, nrun);
                        }
                        res1.store_long();
		if(nobjecev > 100){
			res1.store();
		}

                        if(iprint >= 3) {
                            cout  << "  NEW OPTIMUM" << endl;
                               prt8(n, vm, xopt, x, des, report_all, direct);
                        }
                        xopt = xp;
                        fopt = fp;
                        momopt = momv;
                        nnew++;
                    }
                }
                else{
                //  if the point is lower, use the Metropolis criteria to decide on
                //  acceptance or rejection.
                    p = exp((fp - f)/t);
		if(nobjecev > 100){
			res1.store();
		}
                    pp = ranmar();
                    if (pp < p) {
                        if(iprint >= 3)   prt6(max);
                        x = xp;
                        momd = momv;
                        f = fp;
                        nacc = nacc++;
                        nacp(h) = nacp(h) + 1;
                        ndown++;
                    }
                    else{
                        nrej++;
                        if(iprint >= 3)   prt7(max);
                    }
                }
            }
        }

	cout << "Adjust VM" << endl;
        //  Adjust VM so that approximately half of all evaluations are accepted.
        for(int i = 0; i < nopt; i++){
            ratio = ((double) nacp(i)) / ((double) ns);
            if (ratio > .6) {
                vm(i) = vm(i)*(1.0 + c(i)*(ratio - 0.6)/0.4);
            }
            else if (ratio < 0.4) {
                    vm(i) = vm(i)/(1.0 + c(i)*((0.4 - ratio)/0.4));
            }
            if (vm(i) > (ub(i)-lb(i))) {
                    vm(i) = ub(i) - lb(i);
            }
        }
        if(iprint >= 2) {
           prt8(n, vm, xopt, x, des, report_all, direct);
        }
        
        for(int j = 0; j < nopt; j++){
            nacp(j) = 0;
        }
        
    }
    if(iprint >= -1000000000) {
        prt9(max,n,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew,
             nobjecev,momopt,nmom, des, direct, res1.vimmortal, res1.ns);
    }
    
    //  Check termination criteria.
    quit = false;
    fstar(0) = f;
    diff = fopt - fstar(0);
    if (abs(diff) <= eps) {
        quit = true;
    }
    for(int i = 1; i < neps; i++){
        diff = f - fstar(i);
        if (abs(diff) > eps) quit = false;
    }

    //  Terminate SA if appropriate.
    if (quit) {
        x = xopt;
        ier = 0;
        if (! max) fopt = -fopt;
        if(iprint >= 1)   prt10();
        goto end;
      
    }
    
    //  if termination criteria is not met, prepare for another loop.
    t *= rt;
    for(int i = neps-1; i > 0; i--){
        fstar(i) = fstar(i-1);
    }
    f = fopt;
    x = xopt;
    momd = momopt;
    goto max_loop;
    
    end:
        res1.wipe_long();
        cout << "LOWEST VALUE: " << fopt << endl;
        momv = momopt;
        fret = fopt;
        return xopt;
        
}
