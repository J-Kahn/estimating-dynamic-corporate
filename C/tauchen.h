#include <cmath>
#include "matrix.h"

#ifndef TAU_H
#define TAU_H

using namespace std;

double normalCDF(double value);

double normalPDF(double x, double mu, double s2);
 
void shussey(double theta, double rho, double sigma, matrix y, matrix &A, matrix &P);

void tauchen(double theta, double rho, double sigma, double sigmaz, int n, double* A, double* P);

void thussey(double theta, double rho, double sigma, double sigmaz, int n, matrix &A, matrix &P);
#endif