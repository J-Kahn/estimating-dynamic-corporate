#include <cmath>
#include <omp.h>
#include "tauchen.h"

using namespace std;


#ifndef SOLVE_H
#define SOLVE_H

#define _USE_MATH_DEFINES





void intPrint(int *p, string file, int rows, int cols = 0);

class dynprob{
    public:

        // Vector of parameters
        double *pars, *scale;
        int *porder;
        int thr_nums;
        // Policy functions
        double* vold, *vnew, *vancient, *vimmortal, *scp, *rp, *kp, *bp, *dp, *mktbp, *As,*bs,* ks, *divs, *pay;


        // Transistion matrix and grid for Z
        double* P, *A;
        
        // Dividends, expectation of V, payoff grid
        double *scpart, *Evold;// *pay;

        // Grids for state variables
        double* b, *k, *cp;
        double kss;
        // Grids for policy variables
        double *bint, *kint, *nintover;

        // Interpolate weights
        double *wbint, *wkint;

        // Prices
        double *q, *wpart;

        // Policy indexes
        int *pnew, *iindex, *bindex;
        
        double tauch_sds;
        int delp, s_int_per_thread_num;
        int storer, override_howard, override_mcpd;

        int iter; // Keep track of stuff during solution
        int nk, nb, nbint; // Sizes of the problem
        int  nc, na, ns, nkint, nt, nper, bzero; // A bunch of other stuff
        double r_f; // interest rates
        double phi, xhi, xlo;
	    double rho, mu, sigma, beta, delta, alpha, lambda0, lambda1;
	    double cost_k_up, cost_k_down, smoothc_up;
	    double lambda2, flowcost, r_p, firesale, fixc, smoothc, tau_c;
	    double debtc_up, debtc_down, debtc_fixed_up, debtc_fixed_down;
	    double debtc_adj_up, debtc_adj_down, tau_d, cashc;
        double lambda_up, lambda_d1, lambda_d2; // issuance costs/dividend taxes
 
        // Initialize model class by the sizes of state space grid
        dynprob(int nbl, int nbpl, int nkl, int nkpl, int nal);
        dynprob(int nbl, int nbpl, int nkl, int nal);
        
        void grid_setup();

        // Just a wrapper for the tauchen function in tauchen.cpp
	    void lom_tauchen(double muz, double rhoz, double sigmaz){

	        mu = muz, rho = rhoz, sigma = sigmaz;
	         
	        tauchen(mu, rho, sigma, sigma * tauch_sds, na, A, P);
	    }

	    void lom_tauchen(){
	        
	        mu = pars[21], rho = pars[22], sigma = pars[23];
	        tauchen(mu, rho, sigma, sigma * tauch_sds, na, A, P);
	    }
	    double int_interpolate(int *kpts, int b_low, double wb, int i_p, double wz);
	    
	    	    //void int_interpolate(int *kpts, int b_low, double wb, int i_p, double wz);
      
        /*-----------------------------------------------------------
            Generate profit and expected value (this model only)
            Takes as given a wage grid and risky interest rate grid
            Wage grid and interest rate grid should be nc * ns
        ------------------------------------------------------------*/
        void flow();

        // Set up the interpolate grids
        void interpolant();

        double interpolate(double* kpts, int ilow, double wb, double wk);

        void wipe_expected_value(int this_thread, int ca_per_thread, int thr_num){
            int start_val = ca_per_thread * this_thread, end_val = ca_per_thread * (this_thread + 1);
 
            if(this_thread == thr_num - 1){
                end_val = na * nc;
            }

            for(int s = start_val; s < end_val; s++){
                Evold[s] = 0;
            }
        }
        
       // Wipe full payoff (current dividend + future value) by setting it to current dividend
        void wipe_period_value(int this_thread, int s_per_thread){
            int start_val = s_per_thread * this_thread * nc, end_val = nc * s_per_thread * (this_thread + 1);

            for(int s = start_val; s < end_val; s++){
                pay[s] = divs[s];
            }
        }
        
        
        double interpolate(double* kpts, int b_low, double wb, int i_p, double wz);
        double interpolate(double* kpts, int b_low, double wb, int i_p);

        void expected_value(int this_thread, int s_int_per_thread, int thr_num);

        void expected_single_value(double &ev, int ai, int cp);


        void choice(int &c, double &u, int s);

        void howard_update(double& u, int s, int c);

        /*-----------------------------------------------------------
            Solve model
        ------------------------------------------------------------*/
        void vfi(const double tol, const double maxtol, const int maxint, const int verbose, const int McQueen, const int hmiter, const int thr_num, int store = 0){
            // Turn pars into doubles.

            double maxd = 50, hdiff =50, mind = 50000, mcpd = 50, mcpdadd =0, convd = 50;
            int hiter = 0;
            thr_nums = thr_num;
            omp_set_num_threads(thr_num);
            double maxds[thr_num], minds[thr_num];
            
            //matrix vancient = vold;
            
            double total_expected_value_time = 0, total_maximization_time = 0, total_adjustment_time = 0;
            double t = omp_get_wtime(), bigt = omp_get_wtime();
            
            //doubleCopy(vold, vnew, ns);
            int first;
            iter = 0;
            
            if(storer == 0){
              reset();
            }
            
	    wipe();
            int s_per_thread = ns / thr_num;
            int ca_per_thread = na * nc / thr_num ;
            int s_int_per_thread = s_int_per_thread_num / thr_num;

            int this_thread;


            int nh = (int) (5 / (1 - beta));
            int hburn = 10;
            if(override_howard == 0){
              if(hmiter == 0){
                  hburn = 100000;
              }
            }

            #pragma omp parallel private(this_thread)
            {
                this_thread = omp_get_thread_num();

                wipe_expected_value(this_thread, ca_per_thread, thr_num);
		#pragma omp barrier

                expected_value(this_thread, s_int_per_thread, thr_num);
            }



            while(convd>tol){

                #pragma omp barrier


                #pragma omp master
                {
                    if(verbose > 1){
                      total_expected_value_time += omp_get_wtime() - t;
                      cout << "EV time: " << omp_get_wtime() - t << std::endl;
                    }
                    t = omp_get_wtime();
                    
                    for(int i = 0; i < thr_num; i++){
                                     maxds[i] = -INFINITY;
                minds[i] = INFINITY;
                    }
                }
                #pragma omp barrier

                #pragma omp parallel for private(this_thread)
                for(int s = 0; s < ns; s++){
                this_thread = omp_get_thread_num();
                    double v = vold[s];
                    double u = -5000;
                    int c = 0;

                    //payoff_add(s);

                    choice(c, u, s);

                    vold[s] = u;
                    pnew[s] = c;
                                
                    // Update bounds
                    if(u - v > maxds[this_thread]){
                        maxds[this_thread] = u - v;
                    }
                    if(u - v < minds[this_thread]){
                        minds[this_thread] = u - v;
                    }
                }
                
                #pragma omp barrier
                
                #pragma omp master
                {
                    if(verbose>1){
                      total_maximization_time += omp_get_wtime() - t;
                      cout << "Maximization time: " << omp_get_wtime() - t << std::endl;
                    }

                    t = omp_get_wtime();
                    
                    maxd = maxds[0];
                    mind = minds[0];
                    for(int i = 1; i < thr_num; i++){
                      if(maxds[i] > maxd){
                        maxd = maxds[i];
                      }
                      if(minds[i] < mind){
                        mind = minds[i];
                      }
                    }

                    /*tt = clock() - t;
                    if(verbose > 0){
                        Rprintf ("Time to 2 is %f seconds.\n",((float)tt)/CLOCKS_PER_SEC);
                    }*/
                    convd = abs(maxd);
                    if(abs(mind)>maxd){
                        convd = abs(mind);
                    }
                    
                    mcpd = beta/(1-beta)*(maxd-mind);
                    
                    // Update McQueen-Porteus maximum
                    //mcpd = beta / (1 - beta) * (maxd - mind);
                    mcpdadd = beta/(1-beta) * (maxd+mind)/2.0;
                    if(McQueen == 2 || override_mcpd == 2){
                        if(abs(mcpd) < convd){
                           convd = abs(mcpd);
                        }
                        if(iter > 0){
                          doubleAdd(vold, beta/(1-beta) * (maxd+mind)/2.0, ns);
                        }
                    }
                    else if(McQueen == 1 || override_mcpd == 1){
                      if(abs(mcpd) < convd){
                        convd = abs(mcpd);
                      }
                    }
                    
                    if(convd>maxtol){
                          iter = maxint + 10;
                    }

                    iter++;
                    
                    if(iter >= maxint){
                        convd = 0;
                          iter = maxint + 10;
                        //vold = vancient;
                        //vnew = vancient;
                    }

                    if(verbose > 0){
                        cout << "Iteration: " << iter << " of " << maxint << " maxd " << maxd << " mind " << mind << " mcpd " << mcpd << " mcpd_add " << mcpdadd << " convd " << convd << std::endl;
                    }

                    hdiff = 50;
                    hiter = 0;

                    if(verbose>1){
                      total_adjustment_time += omp_get_wtime() - t;
                      cout << "Adjustment time: " << omp_get_wtime() - t << std::endl;
                    }
                }

                if(iter > hburn){
                    for(int hi = 0; hi < nh; hi++){


                        #pragma omp parallel for
                        for(int s = 0; s < ns; s++){
                            double u = 0;
                            int cip = pnew[s];
                            //payoff_add(s);

                            howard_update(u, s, cip);

                            vnew[s] = u;
                        }

                        #pragma omp barrier

                        #pragma omp parallel for
                        for(int s = 0; s < ns; s++){

                            vold[s] = vnew[s];
                        }


                    }
                }

                #pragma omp barrier

                #pragma omp parallel private(this_thread)
                {
                    this_thread = omp_get_thread_num();

                    wipe_expected_value(this_thread, ca_per_thread, thr_num);
		    #pragma omp barrier
                    expected_value(this_thread, s_int_per_thread, thr_num);
                }
                             

            }
            
            #pragma omp barrier

            if(McQueen == 1 || override_mcpd == 1){
                if(mcpd < tol && iter < maxint){
                    doubleAdd(vold, beta/(1-beta) * (maxd+mind)/2, ns);
                }
            }
            if(verbose>=0){
              cout <<  " iter: " << iter;
              cout << " total time: " << omp_get_wtime() - bigt;
              cout << " maximization time: " << total_maximization_time;
              cout << " expectation time: " << total_expected_value_time;
              cout << " adjustment time: " << total_adjustment_time;
              cout << std::endl;
            }
        }
       
        double interpolate(double* kpts, int ilow, double wb, double wk, double wz);     
        void policy();

        void wipe(){
            for(int s = 0; s < ns; s++){
                vold[s] = vancient[s];
            }
        }
        
      void reset(){
            for(int s = 0; s < ns; s++){
                vold[s] = 0.0;
	            	vancient[s] = 0.0;
            }
        }
       void reset_long(){
            for(int s = 0; s < ns; s++){
                vold[s] = 0.0;
	            	vancient[s] = 0.0;
                vimmortal[s] = 0.0;
            }
        }   
        void store(){
            for(int s = 0; s < ns; s++){
                vancient[s] = vold[s];
            }
        }
        // Print out a bunch of csv files for policy functions
        void store_long(){
            for(int s = 0; s < ns; s++){
                vimmortal[s] = vancient[s];
            }
        }
        void wipe_long(){
            for(int s = 0; s < ns; s++){
                vancient[s] = vimmortal[s];
                vold[s] = vimmortal[s];
            }
        }
        void wipe_long(matrix mat){
            for(int s = 0; s < ns; s++){
                vancient[s] = mat(s);
                vimmortal[s] = mat(s);
                vold[s] = mat(s);
            }
        }
        void print();
        void print(string file);
        
        void print_pi(string file);

        void print_vecs(string title);

};

#endif
