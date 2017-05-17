#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "matrix.h"
#include "tauchen.h"
using namespace std;

double normalCDF(double value){
    return 0.5 * erfc(-value * M_SQRT1_2);
}

double normalPDF(double x, double mu, double s2){
    return 1/sqrt(2*3.1415926*pow(s2,2)) * exp(-pow((x-mu)/s2,2)/2);
}
 
void ghquad(int n, int maxit, matrix &x, matrix &w){
/* Adapted from: Press, et al. Numerical Recipes in Fortran, 2nd ed.
 Returns the abscissas (x) and weights (w) for Gauss-Hermite
 quadrature of order n.
 This is used to compute integrals of the form
    int(exp(-x**2)f(x))
 over the whole real line.
 The integral is computed as w'f(x)
 The integration will be exact for any polynomial f(x) of order
 less than 2n. maxit - maximum iteratioins.
 Example:
   To compute the integral of the function exp(-x**2)(x**4+10) use:
   [x,w]=ghquad(3)
   int=w'(x**4+10)

 To convert these abscissas and weights to those appropriate for
 computing moments of the standard normal distribution
 multiply the abscissas by sqrt(2) and divide the weights by
 sqrt(pi). */

    double tol, pim4, z, z1, p1, p2, p3, pp, pdist;
    matrix xtemp(n);
    int m, i, j, its;


    tol=100;
    pim4= 0.7511255444649425;
    m= (n+1)/2;
    w=x;


    for(i = 0; i < m; i++){


        if(i==0){       
            z=sqrt(2.0*((double) n)+1.0-1.85575*pow(2.0*((double) n)+1.0,-1.0/6.0));
        }
        else if(i==1){
           z=z-1.14*pow(((double) n),0.426)/z;
        }
        else if(i==2){
           z=1.86*z-0.86*x(0,0);
        }
        else if(i==3){
           z=1.91*z-0.91*x(1,0);
        }
        else{
           z=2.0*z-x(i-2,0);
        }


        its = 0;
        pdist = 80;
        while(pdist>1e-14 && its < maxit){
           its++;
           p1=pim4;
           p2=0.0;
           for(j=0; j < n; j++){
                   p3=p2;
                   p2=p1;
                   p1=z*sqrt(2.0/((double) j+1))*p2-sqrt((((double) j))/((double) j+1))*p3;
           }
           pp=sqrt(2.0*((double) n))*p2;
           z1=z;
           z=z1-p1/pp;
           pdist = abs(z-z1);


        }
        if(its>=maxit){
            cout << "Exceeded allowable iterations in finding root" << endl;
        }

        x(n-i-1,0)=-z;
        x(i,0)=z;

        w(i,0)=2/(pp*pp);
        w(n-i-1,0)=w(i,0);
    }
}

void shussey(double theta, double rho, double sigma, matrix y, matrix &A, matrix &P){

    // Approximates the definite integral of f from a to b by the composite Simpson's rule, using n subintervals
    int n = y.nrows();
    matrix t1(n);
    matrix t(n,n), t2(n,n);
    double h = y(1) - y(0);
    int iteration, i, j, ii, jj;
                                       // Calculate Gauss Hermite grid and weight
    iteration=20000;                     // Maximum iterations


                                       // Calculate transition matrix
    double exz;
    double w;
    for(i=0; i < n; i++){
        for(j = 0; j < n;  j++){
            exz = (1-rho) * theta + rho * y(i);
            if(j == 0){
                w = 1 * h / 3;
            }
            else if(j == n - 1){
                w = 1 * h / 3; 
            }
            else{
                if(j % 2 == 0){
                    w = 2 * h / 3;
                }
                else{
                    w = 4 * h / 3;
                }
            }
            t(i,j)= normalPDF(y(j), exz, sigma) / w;
        }
    }

    //Normalize each row
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            t1(i) += t(i, j);
        }
    }
    for(ii = 0; ii < n; ii++){
        for(jj=0; jj < n; jj++){
            t2(jj,ii) = t(ii,jj)/t1(ii);
        }
    }
    
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            t(i, j) = t2(j, i);
        }
    }
    
    A = y;
    P = t;
    for(int i = 0; i < n; i++){
        A(i) = exp(y(i)) / (1+exp(y(i)));
    }
    
 }


// Tauchen approximation to an AR 1 process
void tauchen(double theta, double rho, double sigma, double sigmaz, int n, double* A, double* P){
    
      // Setup the Z matrix
    double   m=sigmaz/sigma;

    // Generate the first state
    A[n-1] = m*(sigma/std::sqrt(1-rho*rho));
    A[0] = -m*(sigma/std::sqrt(1-rho*rho));
    double df = (A[n-1]-A[0])/((double) n-1);
    
    // Generate n 2:n
    for(int i = 1; i < n-1; i++){
            A[i] = A[i-1]+df;
    }
    
    //Setup the transition matrix
   int ina;
      
    //Start generating the n, looping over rows then columns.
    for(int i = 0; i < n; i ++){
        ina = i * n;
        for(int j = 0; j < n; j ++){
            if(j==0){
                P[ina + j] = normalCDF((A[0]+df/2-rho*A[i])/sigma);

            }
            else{
                if(j==n-1){
                    P[ina + j] = (1-normalCDF((A[n-1]-df/2-rho*A[i])/sigma));

                } 
                else{
                    P[ina + j] = normalCDF((A[j]+df/2-rho*A[i])/sigma)-normalCDF((A[j]-df/2-rho*A[i])/sigma);
                }
            }
        }
    }
        
    for(int i = 0; i < n; i ++){
        A[i] = std::exp(A[i] + theta);
    }
}

 void thussey(double theta, double rho, double sigma, double sigmaz, int n, matrix &A, matrix &P){
    matrix  y(n), t(n,n), w(n);
    
    int iteration, i, j, ii, jj;
    
    matrix t1(n), x(n), t2(n,n), xtemp(n,1);
 
                                       // Calculate Gauss Hermite grid and weight
    iteration=20000;                     // Maximum iterations
 
    ghquad(n,iteration,x,w);
    y=sqrt(2.0)*sigmaz*x+theta; //(1.0-rho); // Calculate corresponding y
    w = w/sqrt(3.1415926);
                                       // Calculate transition matrix
    double exz;
    for(i=0; i < n; i++){
        for(j = 0; j < n;  j++){
            exz = (1-rho) * theta + rho * y(i);
            t(i,j)=w(j) * normalPDF(y(j), exz, sigma) / normalPDF(y(j),theta, sigmaz);
        }
    }

    //Normalize each row
    t1 = rsum(t);
 
    for(int ii = 0; ii < n; ii++){
        for(int jj=0; jj < n; jj++){
            t2(jj,ii) = t(ii,jj)/t1(ii);
        }
    }
 
    t=t2.t();

    A = exp(y);
    P = t;
    for(int i = 0; i < n; i++){
        A(n - 1 - i) = exp(y(i));
        for(int j = 0; j < n; j++){
            P(n - 1 - i, n - 1 - j) = t(i, j);
        }
    }

 }