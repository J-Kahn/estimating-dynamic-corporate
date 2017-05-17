#include <cmath>
#include <omp.h>
#include "tauchen.h"
#include "solve.h"

using namespace std;

        // Initialize model class by the sizes of state space grid
        dynprob::dynprob(int nbl, int nbpl, int nkl, int nkpl, int nal){
            // What did we just read in?
            nb   = nbl,
            nbint = nbpl,
            na    = nal,
            nk    = nkl,
            nkint = nkpl;

            // What size will the value function, policy function and payoff grid be?
            ns    = nb * nk * na,
            nc    = nkint * nbint;
            tauch_sds = 2;
            delp = 4;
            s_int_per_thread_num = na * nc;

            // Make all the grids the right size
            b      = new double[nb],
            wbint  = new double[nc],
            wkint  = new double[nc],
            bint   = new double[nbint],
            k      = new double[nk],
            kint   = new double[nkint],
            vold   = new double[ns],
            vancient   = new double[ns],
            vimmortal  = new double[ns],
            vnew   = new double[ns],
            q      = new double[nc * na],
            kp     = new double[ns],
            bp     = new double[ns],
            dp     = new double[ns],
            mktbp  = new double[ns],
            scp    = new double[ns],
            rp     = new double[ns],
            wpart  = new double[ns],
            scpart = new double[nk * na * nc],
            Evold  = new double[na * nc]; // note: EXPECTATION OVER Z ONLY;
            P      = new double[na * na];
            A      = new double[na];
            override_howard = 1;
            override_mcpd = 0;
            storer = 1;
            
            // This is now too large.
            pars  = new double[24];
            scale = new double[24];
            porder = new int[24];
            
            // Initialize the parameter vector
            doubleZeros(pars, 24);
            doubleOnes(scale, 24);
            intSeq(porder, 0, 24);
            doubleZeros(vold, ns);
            scale[11] = 10.0;

      			beta             = pars[0],  // beta rate
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
            
      		  r_f = ((1 / beta) - 1)/(1 + tau_c); // Set beta rate
      			
            // And now we need some other things the right size
            iindex  = new int[nc];
            pnew    = new int[ns];
        }

         void dynprob::grid_setup(){
      			beta             = pars[0],  // beta rate
      			alpha            = pars[1],  // curvature
      			delta            = pars[2],  // depreciation
      			tau_c            = pars[14]; // downward debt issuance cost
        			firesale         = pars[9],  // firesale value of capital
    			
          r_f = ((1 / beta) - 1)/(1 + tau_c); // Set beta rate
          double kstar= std::pow((r_f + delta) / (alpha), 1 / (alpha - 1));
          double kmax = std::pow((r_f + delta) / (alpha * A[na - 1]), 1 / (alpha - 1));
          //kmax = kmax * 0.33;
    
          double kmin = std::pow((r_f + delta) / (alpha * A[0]), 1 / (alpha - 1));
          //cout << A[0] << " " << A[na-1] << endl;
          //kmin = kmin * 3;
          
          /*if(kmin > kmax){
            kmax = kmax / 0.33;
            kmin = kmin /3.0;
          }*/

          kmax = std::log(kmax);
          kmin = std::log(kmin);

          for(int i = 0; i < nk; i++){
              k[i] = std::exp(kmin + (kmax - kmin) * (double) i / (nk - 1));
            //cout << "k " << k[i] << endl;
              
          }
          
         for(int i = 0; i < nkint; i++){
              kint[i] = std::exp(kmin + (kmax - kmin) * (double) i / (nkint - 1));
          }
    
          double bmax = firesale*exp(kmax);
          double bmin = -bmax;
    
          for(int i = 0; i < nb; i++){
              b[i] = bmin + (bmax - bmin) * (double) i / (nb - 1);
          }
          
         for(int i = 0; i < nbint; i++){
              bint[i] = bmin + (bmax - bmin) * (double) i / (nbint - 1);
          }
        };

        /*-----------------------------------------------------------
            Generate profit and expected value (this model only)
            Takes as given a wage grid and risky interest rate grid
            Wage grid and interest rate grid should be nc * ns
        ------------------------------------------------------------*/
         void dynprob::flow(){
            
			beta             = pars[0],  // beta rate
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
      
			double r_f = ((1 / beta) - 1)/(1 + tau_c); // Set beta rate
			
            #pragma omp parallel for
            for(int ki = 0; ki < nk; ki ++){
                double kbpart =  pow(k[ki],alpha); // Pow function is expensive
                
                for(int ai = 0; ai < na; ai++){
  
                    double pipart = A[ai] * kbpart; // And I'm lazy
                
                    for(int bi = 0; bi < nb; bi++){
                        int sindex = ki * nb * na + bi * na + ai;
                        wpart[sindex] = pipart - (1 + r_f) * b[bi]; // And I don't want to access memory too much
                        //double tpart = pipart;
                        //wpart[sindex] += - tau_c * tpart;
                        /*if(tpart > 0){
                            wpart[sindex] += - tau_c * tpart;
                        }*/
                    }
                }
            }

            #pragma omp parallel for
            for(int ki = 0; ki < nk; ki ++){
                 int sindexc = ki * nc * na;
                 for(int kpi = 0; kpi < nkint; kpi++){
                    for(int bpi = 0; bpi < nbint; bpi++){
                        for(int ai = 0; ai < na; ai++){
                            double cf =  0;
                            int cp = kpi * nbint + bpi;
                            
                            // Investment costs
                            if(k[ki] * (1 - delta) < kint[kpi]){
                                cf += (1 + cost_k_down) * (kint[kpi] - (1 - delta) * k[ki]);
                            }
                            else{
                                cf += (1 + cost_k_up) * (kint[kpi] - (1 - delta) * k[ki]);
                            }
                            
                            if(bint[bpi] > firesale * kint[kpi] * (1 - delta)){
                              cf += 100000 * kint[nkint - 1];
                            }
                            cf += 0.5 * smoothc * (kint[kpi]/k[ki] - (1 - delta)) * (kint[kpi]/k[ki] - (1 - delta));
                                
                            scpart[sindexc + ai * nc + cp] = cf - bint[bpi]; // Store it in a big fat vector
                        }
                    }
                }
            }
        }

        // Set up the interpolate grids
         void dynprob::interpolant(){
            int underk = 0, underb = 0;
            double wbint_temp[nbint], wkint_temp[nkint];
            int bindex[nbint], kindex[nkint];
            string title = "test/test";
            for(int bip = 0; bip < nbint; bip++){
                if(bip < nbint - 1){
                    if(bint[bip] >= b[underb+1]){
                        underb++;
                    }
                }

                bindex[bip] = underb;
                wbint_temp[bip]  = (bint[bip] - b[underb]) / (b[underb+1] - b[underb]);
            }
            
            for(int kip = 0; kip < nkint; kip++){
                if(kip < nkint - 1){
                    if(kint[kip] >= k[underk+1]){
                        underk++;
                    }
                }
                
                kindex[kip] = underk;
                wkint_temp[kip]  = (log(kint[kip]) - log(k[underk])) / (log(k[underk+1]) - log(k[underk]));
            }
            
            
            #pragma omp parallel for
            for(int ki = 0; ki < nkint; ki++){
                for(int bi = 0; bi < nbint; bi++){
                    int si = ki * nbint + bi;
                    iindex[si] = (kindex[ki] * nb + bindex[bi])*na;
                    wbint[si] = wbint_temp[bi];
                    wkint[si] = wkint_temp[ki];
                   // cout << si << " " << iindex[si] << " " << kindex[ki] << " " << nindex[ni] << " " << bindex[bi] << " " << ki << " " << ni << " " << bi << " " << bint[bi] << " " << b[nb - 1] << " " << nbint << " " << nb << endl;
                }
                
            }

        }

         double dynprob::interpolate(double* kpts, int ilow, double wb, double wk){
          
          double k_p = 0;
          //cout << ilow + ((nn + 1) * nb + 1) * na << endl;
          k_p += (1 - wb) * (1 - wk) * (kpts[ilow]);
          k_p += (1 - wb) * wk       * (kpts[ilow + nb * na]);
          k_p += wb       * (1 - wk) * (kpts[ilow + na]);
          k_p += wb       * wk       * (kpts[ilow + (nb + 1) * na]);

          return k_p;
        }

         double dynprob::interpolate(double* kpts, int ilow, double wb, double wk, double wz){
          
          double k_p = 0;
          //cout << ilow + ((nn + 1) * nb + 1) * na << endl;
          k_p += (1 - wb) * (1 - wk) * (kpts[ilow] * (1 - wz) + kpts[ilow + 1] * wz);
          k_p += (1 - wb) * wk       * (kpts[ilow + nb * na] * (1 - wz) + kpts[ilow + nb * na + 1] * wz);
          k_p += wb       * (1 - wk) * (kpts[ilow + na] * (1 - wz) + kpts[ilow + na + 1] * wz);
          k_p += wb       * wk       * (kpts[ilow + (nb + 1) * na] * (1 - wz) + kpts[ilow + (nb + 1) * na + 1] * wz);

          return k_p;
        }

         void dynprob::expected_value(int this_thread, int s_int_per_thread, int thr_num){
            // calculate expected value
            double wb, wk, pt, w1;
            int underi = 0, cs;
            // calculate expected value
            
            int bip, kip, rowi = 0, coli = 0, ai, ap, indi, andi;
            int end_val = s_int_per_thread * (this_thread + 1);

            if(this_thread == thr_num - 1){
                end_val = na * nc;
            }
            
            for(int s = s_int_per_thread * this_thread; s < end_val; s++){
                ai = s / nc;
                cs = s - ai * nc;
                andi = ai * na;
                
                underi = iindex[cs];
                wb     = wbint[cs];
                wk     = wkint[cs];
                
                w1 = 0;
                //cout << underi << " of " << ns <<endl;
                for(ap = 0; ap < na; ap++){
                    w1     += beta * P[andi + ap] * interpolate(vold, underi + ap, wb, wk);
                }
                //cout << "donezo" << endl;
                Evold[s] = w1;
            }
        }

         void dynprob::expected_single_value(double &ev, int ai, int cp){

            int underi = iindex[cp];
            double wb     = wbint[cp];
            double wk     = wkint[cp];
            int andi = ai * na;

            for(int ap = 0; ap < na; ap++){

                ev += beta * P[andi + ap] * interpolate(vold, underi + ap, wb, wk);
            }
        }


         void dynprob::choice(int &c, double &u, int s){

            int ai = s % na;
            int aindex = ai * nc;
            int indi = s * nc;
            int ki   = s / na / nb;
            int si   = (ki * na + ai) * nc;
            int bpi, kci;

            double wi = wpart[s];

            double pi =  wi - scpart[si];

            if(pi < 0){
                pi += - lambda0 + lambda1 * pi;
            }
            else{
                pi += -lambda2 * pi;
            }

            u = pi + Evold[aindex];
            c = 0;
            double up = 0, ev = 0;

            for(int cp = 1; cp < nc; cp++){

                bpi = cp % nbint;
                kci = cp / nbint;

                pi =  wi - scpart[si + cp];

                if(pi < 0){
                    pi += - lambda0 + lambda1 * pi;
                }
                else{
                    pi += - lambda2 * pi;
                }

                up = pi + Evold[aindex + cp];

                if(up > u){
                        c = cp;
                        u = up;
                }
            }
        }

         void dynprob::howard_update(double& u, int s, int c){
            int ai = s % na;
            int aindex = ai * nc;
            int indi = s * nc;
            int ki   = s / na / nb;
            int si   = (ki * na + ai) * nc;
            int bpi, kci;

            double wi = wpart[s];

            double pi =  wi - scpart[si + c];

            if(pi < 0){
                pi += - lambda0 + lambda1 * pi;
            }
            else{
                pi += -lambda2 * pi;
            }

            double ev = 0;

            expected_single_value(ev, ai, c);

            u = pi + ev;
        }
                
         void dynprob::policy(){
            int p, scindex, ai, ki;
            
            for(int s = 0; s < ns; s++){
                /*if(defas[s] == 1){
                    kp[s] = kint[0];
                    bp[s] = 0;
                    np[s] = nint[0];
                    p     = bzero;
                }
                else{*/
                    p     = pnew[s];

                //}

                ai       = s % na;
                ki       = s / na / nb;
                scindex  = ki * nc * na + ai * nc + p;
                rp[s]    = A[ai] * pow(k[ki],alpha);
                scp[s] = scpart[scindex];
                kp[s]    = kint[p / (nbint)];
                bp[s]    = bint[p % nbint];
                mktbp[s] = vold[s] / k[s];
                dp[s]    = wpart[s] - scpart[scindex];
                
                if(dp[s] < 0){
                    dp[s] += - lambda0 + lambda1 * dp[s];
                }
                else{
                    dp[s] += - lambda2 * dp[s];
                }
            }
        }

         void dynprob::print(string file){
            
            ofstream simfile;
            simfile.open(file.c_str());
          
            simfile << "a,b,k,bp,kp,mktbp,cp,vold,scp,wort,rp" << endl;

            for(int i = 0; i < ns; i++){
              simfile << A[i % na] << "," << b[(i / na) % nb] << ",";
              simfile << k[(i / na / nb)] << ",";
              simfile << bp[i] << "," << kp[i] << "," << mktbp[i] << ",";
              simfile << dp[i] << "," << vold[i] << "," << scp[i] << "," << wpart[i] << "," << rp[i] << endl;
            }
            
            simfile.close();
        }
        
         void dynprob::print_pi(string file){
		ofstream simfile;
		simfile.open(file.c_str());
          for(int s = 0; s < ns; s++){
            int ai = s % na;
            int aindex = ai * nc;
            int indi = s * nc;
            int ki   = s / na / nb;
            int si   = (ki * na + ai) * nc;
            int bpi, kci;

            double wi = wpart[s];
            double pi;
	          int c = 0;
            double up = 0, ev = 0;

            for(int cp =  0; cp <  nc; cp++){

                pi =  wi - scpart[si + cp];

                if(pi < 0){
                    pi += - lambda0 + lambda1 * pi;
                }
                else{
                    pi += - lambda2 * pi;
                }

            simfile << pi << ",    ";

            }
		simfile << endl;
		}
	}

         void dynprob::print_vecs(string title){

            //doublePrint(bs, title + "_bs.csv", ns);
            //doublePrint(As, title + "_as.csv", ns);
            doublePrint(vold, title + "_vold.csv", ns);
            doublePrint(kp, title + "_kp.csv", ns);
            doublePrint(bp, title + "_bp.csv", ns);
            doublePrint(dp, title + "_cp.csv", ns);
            doublePrint(k, title + "_k.csv", nk);
            doublePrint(b, title + "_b.csv", nb);
            doublePrint(bint, title + "_bint.csv", nbint);
            doublePrint(kint, title + "_kint.csv", nkint);
            doublePrint(A, title + "_A.csv", na);
            doublePrint(P, title + "_P.csv", na, na);
        }
