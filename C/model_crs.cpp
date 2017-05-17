#include <cmath>
#include <omp.h>
#include "matrix.h"
#include "tauchen.h"
#include "solve.h"

using namespace std;


/*

This class contains all the value function iteration and policy function objects for simulation.

*/
        dynprob::dynprob(int nbl, int nbpl, int nkl, int nal){
            delp = 4;
            tauch_sds = 2;
            // Update size of things
            nb    = nbl,
            nbint = nbpl,
            na    = nal,
            nk    = nkl;
            ns    = nb * na,
            nc    = nk * nbint;
            s_int_per_thread_num = na * nbint;
            storer = 0;
            // Set grid sizes
            b     = new double[nb];
            bint  = new double[nbint],
            wbint  = new double[nbint];
            k     = new double[nk];
            vold  = new double[ns],
            vnew  = new double[ns],
            vancient = new double[ns],
            vimmortal = new double[ns],
            kp    = new double[ns],
            bp    = new double[ns],
            cp    = new double[ns],
            rp    = new double[ns],
            As    = new double[ns],
            bs    = new double[ns];
            pars  = new double[24];
            scale = new double[24];
            porder = new int[24];
            divs  = new double[ns * nc],
            pay   = new double[ns * nc],
            Evold = new double[na * nc],
            P     = new double[na * na];
            A     = new double[na];
            bindex = new int[nbint];
            pnew = new int[ns];
            override_howard = 0;
            override_mcpd = 0;
            storer = 1;
            
            // Initialize the parameter vector
            doubleZeros(pars, 24);
            doubleOnes(scale, 24);
            intSeq(porder, 0, 24);
            doubleZeros(vold, ns);
            doubleZeros(vancient, ns);
            doubleZeros(vimmortal, ns);

            // Set up parameters
            

      			beta             = pars[0],  // discount rate
      			alpha            = pars[1],  // curvature
      			delta            = pars[2],  // depreciation
      			lambda0          = pars[3],  // fixed issuance cost
      			lambda1          = pars[4],  // lambda
      			lambda2          = pars[5],  // quadratic issuance cost
      			cost_k_down      = pars[6],
      			flowcost         = pars[7],  // flow fixed cost
      			cost_k_up        = pars[8],
      			firesale         = pars[9],  // firesale value of capital
      			fixc             = pars[10], // fixed adjustment cost
      			smoothc          = pars[11], // smooth adjustment cost
      			smoothc_up       = pars[12],
      			tau_d            = pars[13]; // distribution tax rate
      			tau_c            = pars[14]; // downward debt issuance cost
            debtc_up         = pars[15]; // upward debt issuance cost
      			cashc            = pars[16]; // cost for cash balances
      			debtc_fixed_down = pars[17]; // fixed debt adjustment cost
      			debtc_fixed_up   = pars[18]; // fixed debt adjustment cost
      			debtc_adj_up     = pars[19]; // downward debt adjustment cost
      			debtc_adj_down   = pars[20]; // upward debt adjustment cost
            mu               = pars[21];
            rho              = pars[22];
            sigma            = pars[23];

            scale[21] = 10.0;
            scale[11] = 100.0;
            
        }
      
        void dynprob::grid_setup(){
      			firesale         = pars[9];  // firesale value of capital
      			delta            = pars[2];  // depreciation

          double mini = delta * (1.0 - ceil(((double) nk) / (2 * delp))) + delta;
      	
        	for(int i = 0; i <nk; i++){
        		k[i] = mini + ((double) i) * delta / delp;
        		//cout << 'K ' << k[i] << endl;
        	}
        	
          double bhigh = firesale;
      
          for(int i = 0; i < nb; i++){
              b[i] = bhigh * 2.0 * (double) i / (nb - 1) - bhigh;
              //cout << b[i] << endl;
          }
          for(int i = 0; i < nbint; i++){
              bint[i] = bhigh * 2.0 * (double) i / (nbint - 1) - bhigh;
          }
          
        }
      

        // Update A and P from productivity process parameters using Tauchen.
      
      
        /*-----------------------------------------------------------
	        Generate profit and expected value (this model only)
        ------------------------------------------------------------*/

        // Set up flow grid
        void dynprob::flow(){
			beta             = pars[0],  // discount rate
			alpha            = pars[1],  // curvature
			delta            = pars[2],  // depreciation
			lambda0          = pars[3],  // fixed issuance cost
			lambda1          = pars[4],  // lambda
			lambda2          = pars[5],  // quadratic issuance cost
			cost_k_down      = pars[6],
			flowcost         = pars[7],  // flow fixed cost
			cost_k_up        = pars[8],
			firesale         = pars[9],  // firesale value of capital
			fixc             = pars[10], // fixed adjustment cost
			smoothc          = pars[11], // smooth adjustment cost
			smoothc_up       = pars[12],
			tau_d            = pars[13]; // distribution tax rate
			tau_c            = pars[14]; // downward debt issuance cost
      debtc_up         = pars[15]; // upward debt issuance cost
			cashc            = pars[16]; // cost for cash balances
			debtc_fixed_down = pars[17]; // fixed debt adjustment cost
			debtc_fixed_up   = pars[18]; // fixed debt adjustment cost
			debtc_adj_up     = pars[19]; // downward debt adjustment cost
			debtc_adj_down   = pars[20]; // upward debt adjustment cost
      mu               = pars[21];
      rho              = pars[22];
      sigma            = pars[23];
      
			double r_f = ((1 / beta) - 1) / (1 + tau_c); // Set discount rate

            // Parallelize the grid setup
            #pragma omp parallel for
			for(int j = 0; j<nb; j++){
			    for(int m = 0; m < na; m++){
                	double kbpart = A[m]
                   		   - (1 + r_f) * b[j] - flowcost;
                   		   
                    if(b[j] < 0){
                      kbpart += -cashc;
                    }
                    
                    
                   	int sindex = (j * na + m) * nc;
			        for(int l = 0; l < nk; l++){
			            for(int s = 0; s < nbint; s++){

                            // Part of cash flow due to capital / debt choice
				            double cf =  kbpart
                               + bint[s] * (1 - delta + k[l])
                               - 0.5 * smoothc * std::pow(k[l], 2)
                               - k[l] ;


                            // Add fixed cost of capital adjustment
                            if(k[l] < -1e-8){
                                cf += - fixc + cost_k_down * k[l];
                            }
                            if(k[l] > 1e-8){
                                cf += -fixc - cost_k_down * k[l];
                            }
                    
                    if(bint[s] > 0){
                      if(bint[s] > b[j]){
                        if(b[j] < 0){
                          cf += - debtc_fixed_up - debtc_adj_up * (bint[s]) * (bint[s]) - debtc_up * (bint[s]);
                        }
                        else{
                          cf += - debtc_fixed_up - debtc_adj_up * (bint[s] - b[j]) * (bint[s] - b[j]) - debtc_up * (bint[s] - b[j]);
                        }
                      }
                    }
                    
                    
				            // Equity issuance costs
				            if(cf<0){
					            cf += lambda1 * cf - lambda2 * cf * cf / 2 - lambda0;
				            }
				            
				            
				            else{
				                cf *= (1 - tau_d);
				            }
				            
                    
                            // Store in dividends vector

				            divs[sindex + s * nk + l]    = cf;

				        }
			        }
			    }
	        }
        }


        // Set up interpolation for the VFI
        void dynprob::interpolant(){

            // This is just the integer for the lowest bond grid point below the interpolation grid point
            int underb = 0;
            
            // Now we'll go through the (ordered) bond grid,
            for(int bip = 0; bip < nbint; bip++){

                // If you've reached the last grid point, stay below
                if(bip < nbint - 1){

                    // If the next grid point is below the interpolation grid point, move it up
                    if(bint[bip] > b[underb+1]){
                        underb++;
                    }

                }
                
                // store the grid point * na, because that's the index we'll eventually need
                bindex[bip] = underb * na;

                // Store interpolation weight
                wbint[bip]  = (bint[bip] - b[underb]) / (b[underb+1] - b[underb]);
            }
        }

        // THIS IS DEPRECATED
        // Calculate period payoff adding next period expected value to dividend
        /*void dynprob::period_payoff_sum(int this_thread, int s_int_per_thread){
            // calculate expected value
            double wb, w1, w2, pt, fp;
            int underb = 0, underk = 0;
            // calculate expected value
            
            int bip, kip, rowi = 0, coli = 0, ai, ap, indi, andi;
            
            for(int s = s_int_per_thread * this_thread; s < s_int_per_thread * (this_thread + 1); s++){
                ai = s / nbint;
                bip = s - ai * nbint;
                indi = s * nk;
                andi = ai * na;

                w1 = 0;

                underb = bindex[bip];
                wb     = wbint[bip];

                for(ap = 0; ap < na; ap++){
                    w1     += beta * P[andi + ap] * ((1 - wb) * vold[underb+ap] + wb * vold[underb+ap+na]);
                }

                fp = w1 * (1 - delta);
                
                for(kip = 0; kip < nk; kip++){
                    w2 = fp + w1 * k[kip];

                    for(int bi = 0; bi < nb; bi++){
                        pay[bi * na * nc + indi + kip] += w2;
                    }
                }
            }
            
        }


        // THIS IS DEPRECATED
        // Calculate period payoff adding next period expected value to dividend
        void dynprob::period_payoff(int this_thread, int s_int_per_thread){
            // calculate expected value
            double wb, w1, pt, w2;
            int underb = 0, underk = 0;
            // calculate expected value
            
            int bip, kip, rowi = 0, coli = 0, ai, ap, indi, andi;
            
            for(int s = s_int_per_thread * this_thread; s < s_int_per_thread * (this_thread + 1); s++){
                ai = s / nbint;
                bip = s - ai * nbint;
                indi = s * nk;
                andi = ai * na;

                underb = bindex[bip];
                wb     = wbint[bip];
                for(ap = 0; ap < na; ap++){

//                      wkint = (kint(kip) - k(underk)) / (k(underk+1) - k(underk));
                    w1     = beta * P[andi + ap] * ((1 - wb) * vold[underb+ap] + wb * vold[underb+ap+na]);
                    
                    for(kip = 0; kip < nk; kip++){

                        w2 = w1 * (1 - delta + k[kip]);

                        for(int bi = 0; bi < nb; bi++){
                            pay[bi * na * nc + indi + kip] += w2;
                        }
                    }                      //cout << kip << std::endl;
                }
            }
            
        }

        // THIS IS DEPRECATED
        // Maximizes and calculated expectations at the same time
        void dynprob::max_payoff(int this_thread, int s_int_per_thread){
            // calculate expected value
            double wb, w1, pt, w2, up, u;
            int underb = 0, underk = 0;
            // calculate expected value
            
            int bip, kip, rowi = 0, coli = 0, ai, ap, indi, andi;
            
            for(int s = s_int_per_thread * this_thread; s < s_int_per_thread * (this_thread + 1); s++){
                ai = s / nbint;
                bip = s - ai * nbint;
                indi = s * nk;
                andi = ai * na;

                underb = bindex[bip];
                wb     = wbint[bip];

                w1 = 0;

                for(ap = 0; ap < na; ap++){

//                      wkint = (kint(kip) - k(underk)) / (k(underk+1) - k(underk));
                    w1     += beta * P[andi + ap] * ((1 - wb) * vold[underb+ap] + wb * vold[underb+ap+na]);
                }

                for(kip = 0; kip < nk; kip++){

                    w2 = w1 * (1 - delta + k[kip]);

                    for(int bi = 0; bi < nb; bi++){

                        up = divs[bi * na + ai] + w2;
                        u  = vnew[s];

                        if(up > u){
                            vnew[s] = up;
                            pnew[s] = bip * nk + kip;
                        }
                    }
                }
            }
            
        }*/


        // Calculate the expected value in parallel
        void dynprob::expected_value(int this_thread, int s_int_per_thread, int thr_num){
            // calculate expected value

            double wb, w1, pt, fp;
            int underb = 0, underk = 0;
            // calculate expected value
            
            int bip, kip, rowi = 0, coli = 0, ai, ap, indi, andi;
            int end_val = s_int_per_thread * (this_thread + 1);
           

            if(this_thread == thr_num - 1){
                end_val = na * nbint;
            }
            
            for(int s = s_int_per_thread * this_thread; s < end_val; s++){
                ai = s / nbint;
                bip = s - ai * nbint;
                indi = s * nk;
                andi = ai * na;

				underb = bindex[bip];
				wb     = wbint[bip];

                w1 = 0;

                for(ap = 0; ap < na; ap++){

					w1     += beta * P[andi + ap] * ((1 - wb) * vold[underb+ap] + wb * vold[underb+ap+na]);
                }

                fp = w1 * (1 - delta);
                for(kip = 0; kip < nk; kip++){
                        Evold[kip + indi] += fp + w1 * k[kip];
                }
            }
        }

        /*// This parallelizes within the expected value
        void dynprob::expected_value_parallel(){
            // calculate expected value
            double wb, w1, pt, fp;
            int underb = 0, underk = 0;
            // calculate expected value
            
            int rowi = 0, coli = 0, indi, andi;

            #pragma omp parallel for private(indi, andi, underb, wb, w1, fp)
            for(int bip = 0; bip < nbint; bip++){

                underb = bindex[bip];
                wb     = wbint[bip];

                for(int ai = 0; ai < na; ai++){
                    indi = ai * nc + bip * nk;
                    andi = ai * na;

                    w1 = 0;

                    for(int ap = 0; ap < na; ap++){
    //                      wkint = (kint(kip) - k(underk)) / (k(underk+1) - k(underk));
                        w1     += beta * P[andi + ap] * ((1 - wb) * vold[underb+ap] + wb * vold[underb+ap+na]);
                    }

                    fp = w1 * (1 - delta);
                    for(int kip = 0; kip < nk; kip++){
                            Evold[kip + indi] += fp + w1 * k[kip];
                    }                      //cout << kip << std::endl;
                }
            }
        }*/


        // Expected value for a single choice and a
        void dynprob::expected_single_value(double &ev, int ai, int cp){

            int bip = cp % nbint;
            int underb = bindex[bip];
            double wb     = wbint[bip];
            int kip = cp / nbint;
            int andi = ai * na;

            for(int ap = 0; ap < na; ap++){

                ev += beta * P[andi + ap] * ((1 - wb) * vold[underb*na+ap] + wb * vold[underb*na+ap+na]);;
            }

            ev *= (1 - delta + k[kip]);
        }

        // Add up dividends and expected future value DEPRECATED
        /*void dynprob::payoff_add(int s){

            int ai = s % na;
            int aindex = ai * nc;
            int indi = s * nc;

            for(int cp = 0; cp < nc; cp++){
                pay[indi + cp]= divs[indi + cp] + Evold[aindex + cp];
            }
        }*/


        // Select best u by c given s, actual maximization part
        void dynprob::choice(int &c, double &u, int s){

            int ai = s % na;
            int aindex = ai * nc;
            int indi = s * nc;
            u = divs[indi] + Evold[aindex];
            c = 0;
            double up = 0, ev = 0;

            for(int cp = 1; cp < nc; cp++){

                up = divs[indi + cp] + Evold[aindex + cp];

                if(up >= u){
                        c = cp;
                        u = up;
                }
            }
        }



        // Interpolate double variable kpts over bond and z using bond lower point b_low, weight wb,
        // z lower point i_p and weight wz
		double dynprob::interpolate(double* kpts, int b_low, double wb, int i_p, double wz){
		  
		  double k_p = 0;

		  k_p += ((wb)*kpts[(b_low+1)*na + i_p]
		       + (1-wb)*kpts[b_low * na + i_p]) * (1 - wz);

		  k_p += ((wb)*kpts[(b_low+1)*na + i_p + 1]
		       + (1-wb)*kpts[b_low * na + i_p + 1]) * (wz);

		  return k_p;
		}
		
		double dynprob::interpolate(double* kpts, int b_low, double wb, int i_p){
		  
		  double k_p = 0;

		  k_p += ((wb)*kpts[(b_low+1)*na + i_p]
		       + (1-wb)*kpts[b_low * na + i_p]);

		  return k_p;
		}
		
         void dynprob::howard_update(double& u, int s, int c){
        }

        // Interpolate integer variable kpts over bond and z using bond lower point b_low, weight wb,
        // z lower point i_p and weight wz
        double dynprob::int_interpolate(int *kpts, int b_low, double wb, int i_p, double wz){
          
          double k_p = 0;

          k_p += ((wb)* (double) kpts[(b_low+1)  *na + i_p]
               + (1-wb) * (double) kpts[b_low * na + i_p]) * (1 - wz);

          k_p += ((wb) * (double) kpts[(b_low+1)*na + i_p + 1]
               + (1-wb) * (double) kpts[b_low * na + i_p + 1]) * (wz);

          return k_p;
        }

        /*int int_index(int b_low, double wb, int i_p, double wz){

          return b_low * na + i_p;
        }*/
		        
        // Set policy functions from state and pnew index
        void dynprob::policy(){
            int p;
            
            for(int s = 0; s < ns; s++){
                
				p     = pnew[s];
				
                // Get policy functions from state index

				kp[s] = k[p % nk];

				bp[s] = bint[p /nk];

                rp[s] = A[s % na];

                // Dividends from policy and state combination

				cp[s] = divs[s * nc + p];

                // State variable at this point in the grid

				bs[s] = b[s / na];

				As[s] = A[s % na];
            }
        }
        
        void dynprob::print(){

            doublePrint(bs, "bs.csv", ns);
            doublePrint(As, "as.csv", ns);
            doublePrint(vold, "vold.csv", ns);
            doublePrint(kp, "kp.csv", ns);
            doublePrint(bp, "bp.csv", ns);
            doublePrint(cp, "cp.csv", ns);
            doublePrint(rp, "rp.csv", ns);

        }

        // Print out a bunch of csv files for policy functions
         void dynprob::print(string file){
            
            ofstream simfile;
            simfile.open(file.c_str());
          
            simfile << "a,b,bp,kp,cp,vold,rp" << endl;

            for(int i = 0; i < ns; i++){
              simfile << A[i % na] << "," << b[i / na] << ",";
              simfile << bp[i] << "," << kp[i] << ",";
              simfile << cp[i] << "," << vold[i] << "," << rp[i] << endl;
            }
            
            simfile.close();
        }


