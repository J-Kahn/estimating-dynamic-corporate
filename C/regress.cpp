#include "matrix.h"
#include "epf.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


matrix tensor_quad(matrix x, matrix y, int degree){
    matrix xyspline(x.nrows(),6);
    #pragma omp parallel for
    for(int k = 0; k < x.nrows(); k++){
        xyspline(k,0) = 1;
        xyspline(k,1) = x(k);
        xyspline(k,2) = y(k);
        xyspline(k,3) = x(k)*y(k);
        xyspline(k,4) = x(k)*x(k);
        xyspline(k,5) = y(k)*y(k);
    }
    return xyspline;
}

matrix tensor_cube(matrix x, matrix y, int degree){
    matrix xyspline(x.nrows(),10);
    #pragma omp parallel for
    for(int k = 0; k < x.nrows(); k++){
        xyspline(k,0) = 1;
        xyspline(k,1) = x(k);
        xyspline(k,2) = y(k);
        xyspline(k,3) = x(k)*y(k);
        xyspline(k,4) = x(k)*x(k);
        xyspline(k,5) = y(k)*y(k);
        xyspline(k,6) = y(k)*y(k)*x(k);
        xyspline(k,7) = x(k)*x(k)*y(k);
        xyspline(k,8) = x(k)*x(k)*x(k);
        xyspline(k,9) = y(k)*y(k)*y(k);

    }
    return xyspline;
}


matrix tiles(matrix x, int n){
    matrix qrtl(n+1);
    imat numers(x.nrows());
    for(int i = 0; i < x.nrows(); i++){
        numers(i) = i;
    }

    matrix xs = quickOrder(x,numers,0,x.nrows()-1).vect;
    int nx = xs.nrows();

    for(int i = 0; i < n; i++){
        qrtl(i) = xs(i * nx / n);
    }

    for(int i = n-2; i >0; i--){
        if(qrtl(i)==qrtl(i-1)){
            for(int j = 0; j < 10; j ++){
            }
        }
    }

    return qrtl;

}


matrix regress(matrix &y, matrix &x_prime){
    matrix xx = x_prime.cross();
    matrix xy = x_prime % y;

    matrix beta = solve_sym(xx, xy);


    return beta;
}

matrix regress_analysis(matrix y, matrix x_prime, matrix beta, matrix &res, matrix &pred, double &R2){
    int n = y.nrows(), p = x_prime.ncols();
    res = (y-x_prime*beta);
    R2 = 1 - var(res) * ((double) n-1) /((double )n-p) / var(y);
    pred = x_prime*beta;
    beta.print();

    return beta;
}

matrix regress_res(matrix &y, matrix &x_prime, matrix &res){
    int n = y.nrows(), p = x_prime.ncols();
    matrix beta = solve_sym(x_prime.cross(),x_prime % y);
    res = (y-x_prime*beta);
    return beta;
}
