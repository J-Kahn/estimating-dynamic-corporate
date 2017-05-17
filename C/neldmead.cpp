#include <cmath>
#include <omp.h>
#include "matrix.h"
#include "epf.h"
#include "opt.h"
#include <iomanip>

#define _USE_MATH_DEFINES
using namespace std;

// Refined objective functions to check boundaries
double objec_nm(matrix par, matrix m, matrix W, vfiout &res1, data &dat1, matrix &randus, 
                matrix minb, matrix maxb, matrix &mn, 
                    double (*objec) (matrix, matrix, matrix, vfiout&, data&, matrix&, int, matrix&)){
    
    int n = par.nrows();
    int breaker = 0;
    
    double objr;
    
    for(int i = 0; i < n; i ++){
        if(par(i) < minb(i)){
            breaker = 1;
        }
        if(par(i) > maxb(i)){
            breaker = 1;
        }
    }
    
    if(breaker == 1){
        objr = INFINITY;
    }
    else{
        objr = objec(par, m, W, res1, dat1, randus, 10, mn);
    }
    
    return objr;
}


// Nelder-Mead

matrix NelderMead(matrix init_params, int maxevals, double tol, matrix mn, matrix W, int seedn, matrix &f,
                      matrix lb, matrix ub, vfiout &res1, matrix &randus, int nob, string des,
                      double (*objec) (matrix, matrix, matrix, vfiout&, data&, matrix&, int, matrix&)){             
    
    
    int nid;
    int np;

    double maxd = 500;
    int nrun = randus.ncols();
    int nsim = randus.nrows();
    
    int n=init_params.nrows();
    int ns=n+1;
    int shrink=0;
    
    data dat;
    dat = dat_init(nsim, nsim - nob);


    matrix momv;
    matrix params=init_params.col(0),y(n),xbar(n);

    imat order(n+1), numers(n+1);
    double rho=1,chi=2,psi=0.9,sigma=0.9;
    double fr,fe;
    
    ord ranks; 
           

    matrix parS(n,n+1),parSc(n,n+1);
    
    matrix beta;
    
    printf("Evaluate initial points...\n\n");
    double miners = 5000000000000;
    
    for(int i=0;i<n+1;i++){
        cout << i << " of " << n+1 << " ";
        
        numers(i)=i;
        y=init_params.col(i);
        y.print();


        f(i) = objec_nm(y, mn, W, res1, dat, randus, lb, ub, momv, objec);
        parS.col_sub(y, i);
        if(f(i)<miners){
            beta = momv;
            miners = f(i);
        }        
    }
    numers(n) = n;
    f.print();
    parS.print();
    printf("Start Nelder-Mead\n\n");
    ranks = quickOrder(f,numers,0,n);
    order = ranks.index;
    f = ranks.vect;
    order.print();
    for(int i=0;i<n+1;i++){
        parSc.col_sub(parS.col(order(i)), i);
    }
    parS=parSc;
    f.print();
    parS.print();    
    int evals=n+1;
    double maxval=5000;
    
    while(maxval>tol){
        
        shrink=0;
        
        for(int i=0;i<n;i++){
            xbar(i)=0;
            for(int j=0;j<n+1;j++){
            numers(i)=i;
                xbar(i) += parS(i,j) / ((double) n+1);
            }
        }
        
        // Reflection point
        params = (1 + rho) * xbar - rho * parS.col(n);

        fr = objec_nm(params, mn, W, res1, dat, randus, lb, ub, momv, objec);
        evals++;
        
        if(fr<f(0)){
            beta = momv;
            // Expansion point
            y=(1-chi)*xbar+chi*params;
            fe = objec_nm(y, mn, W, res1, dat, randus, lb, ub, momv, objec);
            evals++;
            
            if(fe<fr){
                beta=momv;
                printf("\nEXPAND!\n");
                // Expand
                parS.col_sub(y, n);
                f(n)=fe;
            }
            else{
            printf("\nREFLECT!\n");
                // Reflect
                parS.col_sub(params, n);
                f(n)=fr;
            }
        }
        else{
            if(fr<f(n-1)){
                printf("\nREFLECT!\n");
                // Reflect
                parS.col_sub(params, n);
                f(n)=fr;
            }
            else{
                if(fr<f(n)){
                    
                    // Contraction point
                    y=(1-psi)*xbar+psi*params;
                    fe = objec_nm(y, mn, W, res1, dat, randus, lb, ub, momv, objec);
                    
                    evals++;
                    
                    if(fe<=fr){
                        // Contract outside
                        printf("\nCONTRACT OUTSIDE!\n");
                        parS.col_sub(y,n);
                        f(n)=fe;
                    }
                    else{
                        shrink=1;
                    }
                }
                else{
                    
                    // Contraction point
                    y=(1-psi)*xbar+psi*parS.col(n);
                    fe = objec_nm(y, mn, W, res1, dat, randus, lb, ub, momv, objec);
                    
                    evals++;
                    
                    if(fe<f(n)){
                        printf("\nCONTRACT INSIDE!\n");
                        // Contract inside
                        parS.col_sub(y, n);
                        f(n)=fe;
                    }
                    else{
                        shrink=1;
                    }
                }
                if(shrink==1){
                printf("\nSHRINK!\n");
                    for(int j=1;j<n+1;j++){
                            parS.col_sub(parS.col(0)+sigma*(parS.col(j)-parS.col(0)), j);
                            params=parS.col(j);
                            f(j) = objec_nm(params, mn, W, res1, dat, randus, lb, ub, momv, objec);
                        evals++;
                        }
                    }
                }
            }
     f.print();
    parS.print();  
    numers.print();     
            maxval=0;
            ranks = quickOrder(f,numers,0,n);
            order = ranks.index;
            f = ranks.vect;
            order.print();
            for(int i=0;i<n+1;i++){
                parSc.col_sub(parS.col(order(i)), i);
                
                if(abs(f(i)-f(0))>maxval){
                    maxval=abs(f(i)-f(0));
                }
            }
     f.print();
    parS.print();  
    
                maxd = 0;
            for(int i = 0; i < n; i++){
                if(abs(parSc(i,0)-parSc(i,n))>maxd){
                    maxd = abs(parSc(i,0)-parSc(i,n));
                }
            }
            
            if(maxd < tol){
                maxval = 0;
            }          
            if(evals>maxevals){
                maxval=0;
            }
            parS=parSc;
    
        printf("\n---------------------------\n Function evaluations: %i\n Difference: %f\n Current minimum: %f\n Pars:",evals,maxval,f(0));
        ((parSc.col(0)).t()).print();
        (parSc.col(0)).print("pars_" + des + ".csv");

        cout << " Beta: ";
        (beta.t()).print();
        (beta).print("mom_" + des + ".csv");
        cout << " True: ";
        (mn.t()).print();
        printf("---------------------------\n");

    }
    return(parS);
}

