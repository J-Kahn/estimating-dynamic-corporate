#include "matrix.h"
#define _USE_MATH_DEFINES

#ifndef BSPLINE
#define BSPLINE
using namespace std;

matrix regress(matrix &y, matrix &x);

matrix regress_res(matrix &y, matrix &x_prime, matrix &res);

matrix regress_analysis(matrix y, matrix x_prime, matrix beta, matrix &res, matrix &pred, double &R2);

matrix tensor_quad(matrix x, matrix y, int degree);
matrix tensor_cube(matrix x, matrix y, int degree);

#endif
