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
#define _USE_MATH_DEFINES

/*
  Take jacobians, two estimators at a time, using a linear approximation to
  the Taylor expansion.
*/

void jacobians(matrix par, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0, int l){


    double gamma = 2;

    matrix mn, mn2;

    res1.wipe_long();

    out_auxiliaries(par, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, mn, mn2);
    cout << "Done " << endl;


    matrix report_up, report_down, report, report2, report_up2, report_down2;
    report = mn;
    report2 = mn2;

    matrix v = par * dh0;

    matrix pard = par;
    matrix paru = par;

    int nmom  = mn.nrows();
    int nmom2 = mn2.nrows();
    int n    = par.nrows();
    if(l == -1){
        l = n;
    }

    matrix ys(neval * l *2, nmom);
    matrix ys2(neval * l *2, nmom2);
    matrix xs(neval * l * 2, l * order);

    double dh;
    double fact;
    int error, error_lag, which1, subiter, subiter_lim = 4;

    for(int i = 0; i < l; i++){
        error = 0;
        error_lag = 0;

        dh = v(i) / gamma;
        for(int j = 0; j < neval; j++){
            dh = dh * gamma;
            pard = par;
            paru = par;

            pard(i) -= dh;
            paru(i) += dh;
            res1.wipe_long();

            out_auxiliaries(paru, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, report_up, report_up2);


            if(res1.iter >= 3000){
                error++;
                which1 = 0;
            } else {
        	    if(isnan(report_up(0))){
                  error++;
                  which1 = 0;
        	    } else {
          	    if(isnan(report_up2(0))){
                    error++;
                    which1 = 0;
          	    }
        	    }
            }

       res1.wipe_long();

            out_auxiliaries(pard, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, report_down, report_down2);

            if(res1.iter >= 3000){
                error++;
                which1 = 1;
            } else {
        	    if(isnan(report_down(0))){
                  error++;
                  which1 = 1;
        	    } else {
          	    if(isnan(report_down2(0))){
                    error++;
                    which1 = 1;
          	    }
        	    }
            }

      	    error = 0;
            if(error > error_lag){
              subiter = 0;
                  cout << "Error " << endl;

              while(subiter < subiter_lim){

                  subiter++;
                  error_lag = error;

                  dh /= 10;
                  pard = par;
                  paru = par;
                  pard(i) -= dh;
                  paru(i) += dh;
                  res1.wipe_long();

                  out_auxiliaries(paru, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, report_up, report_up2);

                  if(res1.iter >= 3000){
                      error++;
                      which1 = 0;
                  } else {
              	    if(isnan(report_up(0))){
                        error++;
                        which1 = 0;
              	    } else {
                	    if(isnan(report_up2(0))){
                          error++;
                          which1 = 0;
                	    }
              	    }
                  }

                  res1.wipe_long();

            	    out_auxiliaries(pard, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, report_down, report_down2);

                  if(res1.iter >= 3000){
                      error++;
                      which1 = 1;
                  } else {
              	    if(isnan(report_down(0))){
                        error++;
                        which1 = 1;
              	    } else {
                	    if(isnan(report_down2(0))){
                          error++;
                          which1 = 1;
                	    }
              	    }
                  }
              }

            }

            for(int k = 0; k < nmom; k++){
                ys(i * neval*2 + j*2, k) = report_up(k) - report(k);
                ys(i * neval*2 + j*2 + 1, k) = report(k) - report_down(k);
            }

            for(int k = 0; k < nmom2; k++){
                ys2(i * neval*2 + j*2, k) = report_up2(k) - report2(k);
                ys2(i * neval*2 + j*2 + 1, k) = report2(k) - report_down2(k);
            }

            fact = 1.0;

            for(int ii = 0; ii < order; ii++){
                fact *= (double) (ii + 1);
                xs(i * neval*2 + j*2, i  + ii * l) = std::pow(paru(i) - par(i), ii + 1) / fact;
                xs(i * neval*2 + j*2 + 1, i  + ii * l) = std::pow(par(i) - pard(i), ii + 1) / fact;
            }

            pard = par;
            paru = par;

            error_lag = error;
        }
    }
    cout << error << endl;

    matrix allorder = solve_sym(xs.cross(), xs%ys);

    matrix allorder2 = solve_sym(xs.cross(), xs%ys2);

    jac = matrix(nmom, n);
    jac2 = matrix(nmom2, n);

    for(int i = 0; i < nmom; i++){
      for(int j = 0; j < l; j++){
        jac(i, j) = allorder(j, i);
      }
    }

    for(int i = 0; i < nmom2; i++){
      for(int j = 0; j < l; j++){
        jac2(i, j) = allorder2(j, i);
      }
    }
}


/*
  A simpler method of taking jacobians based on taking two discrete steps. Again
  applied to two estimators at a time to speed up process.
*/
void simple_grad(matrix x0, dynprob &res1, data &dat1, matrix &randus,
matrix (*auxiliary_est)(dynprob &, int, int, data &, matrix&, int), matrix (*auxiliary_est2)(dynprob &, int, int, data &, matrix&, int), int nrun, int thr_num, int order, int neval, matrix &jac, matrix &jac2, double dh0, int l){
    int k = x0.nrows();
    if(l == -1){
        l = k;
    }
    matrix mn, mn2;
    out_auxiliaries(x0, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, mn, mn2);



    int n1 = mn.nrows();
    int n2 = mn2.nrows();
    jac = matrix(n1, k);
    jac2 = matrix(n2, k);

    matrix g1(n1,k), g2(n2,k);
    matrix ax0(k), dh(k), xdup(k), xddw(k);
    matrix argup(k,k), argdw(k,k), f0(n1), f11(n1), f12(n2), arghi(k);
    matrix g1_up(n1,k), g1_dw(n1, k), g2_up(n2,k), g2_dw(n2, k);
    imat isbad(k);
    arghi = x0;
    double dax0;

    int error = 0;

    int redos = 0;
    dax0 = 1.0;
    ax0 = mabs(x0);
    dh = dh0*ax0*dax0;

    for(int ii = 0; ii < k; ii++){
      isbad(ii) = 0;
    }

    dloop:

      xdup = x0 + dh;
      xddw = x0 - dh;

      argup = argup * 0.0;
      argdw = argdw * 0.0;

      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < k; jj++){
          argup(jj, ii) = x0(jj);
          argdw(jj, ii) = x0(jj);
        }
      }

      for(int ii = 0; ii < k; ii++){
        argup(ii,ii) = xdup(ii);
        argdw(ii,ii) = xddw(ii);
      }

      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n1; jj++){
                g1(jj, ii) = 0.0;
        }
      }
      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n2; jj++){
                g2(jj, ii) = 0.0;
        }
      }
      for(int ii = 0; ii < l; ii++){
       res1.wipe_long();
       error = 0;
       for(int jj = 0; jj < k; jj++){
        arghi(jj) = argup(jj,ii);
       }
            out_auxiliaries(arghi, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, f11, f12);

            if(res1.iter >= 3000){
                error++;
                cout << "Iter error: " << res1.iter << endl;
            } else {
        	    if(isnan(f11(0))){
                  error++;
                                    cout << "nAN error 1 " << endl;

        	    } else {
          	    if(isnan(f12(0))){
                    error++;
                    cout << "nAN error 2 " << endl;
          	    }
        	    }
            }

        //make f1 by calling debt with the ith column of argh
      for(int jj = 0; jj < n1; jj++){
        g1_up(jj,ii) = f11(jj);
      }
      for(int jj = 0; jj < n2; jj++){
        g2_up(jj,ii) = f12(jj);
      }
       isbad(ii) = isbad(ii) + error;

        // Repeat with down values
      error = 0;
      for(int jj = 0; jj < k; jj++){
        arghi(jj) = argdw(jj,ii);
      }
         // arghi.print();

       res1.wipe_long();
            out_auxiliaries(arghi, res1, dat1, randus, nrun, thr_num, auxiliary_est, auxiliary_est2, f11, f12);

            if(res1.iter >= 3000){
                error++;
            } else {
        	    if(isnan(f11(0))){
                  error++;
        	    } else {
          	    if(isnan(f12(0))){
                    error++;
          	    }
        	    }
            }


      for(int jj = 0; jj < n1; jj++){
        g1_dw(jj,ii) = f11(jj);
      }
      for(int jj = 0; jj < n2; jj++){
        g2_dw(jj,ii) = f12(jj);
      }
      isbad(ii) = isbad(ii) + error;

    }

    int sumer = 0;
    for(int ii = 0; ii < k; ii++){
      sumer += isbad(ii);
    }
    if(redos < 30){
    if(sumer > 0){
      for(int ii = 0; ii < k; ii++){
        if(isbad(ii) != 0){
          dh(ii) = dh(ii)*0.5;
        }
        isbad(ii) = 0;
      }
      redos++;
      goto dloop;
    }
    else{
      goto end;
    }
    }
    else{
        cout << "FATAL ERROR ESCAPE! " << endl;
        goto end;
    }

    end:
      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n1; jj ++){
          g1(jj, ii) = g1_up(jj, ii) - g1_dw(jj, ii);
        }
      }

      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n2; jj ++){
          g2(jj, ii) = g2_up(jj, ii) - g2_dw(jj, ii);
        }
      }

      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n1; jj ++){
          jac(jj, ii) = g1(jj, ii) / 2.0 / dh(ii);
        }
      }

      for(int ii = 0; ii < k; ii++){
        for(int jj = 0; jj < n2; jj ++){
          jac2(jj, ii) = g2(jj, ii) / 2.0 / dh(ii);
        }
      }
}
