/*
  A bunch of matrix functions, these are wrappers over double pointers
*/


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <omp.h>
#ifndef MAT_H
#define MAT_H

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


typedef std::numeric_limits< double > dbl;


using namespace std;

void doubleCopy(double* matrix_old, double* matrix_new, int ns);

void doubleCopy(double* matrix_old, double* matrix_new, int start, int ns);

void doubleZeros(double* matrix_new, int ns);
void doubleOnes(double* matrix_new, int ns);
void intSeq(int* matrix_new, int ns, int nf);


void doubleZeros(double* matrix_new, int start, int ns);

void doubleAdd(double* matrix_new, double added, int ns);

void doubleAdd(double* matrix_new, double added, int start, int ns);

void doublePrint(double* p, string file, int rows, int cols);

void doublePrint(double* p, string file, int rows);

// Declarations
class matrix;
class imat;
double det(const matrix& a);
matrix diag(const int n);
matrix diag(const matrix& v);
matrix ones(const int rows, const int cols);
int size(const matrix& a, const int i);
matrix zeros(const int rows, const int cols);


double hp_autocorr(matrix mat, matrix matl);
matrix hp_autocorr_infl(matrix mat, matrix matl);
double central_moment(matrix mat, double order);
matrix central_moment_infl(matrix mat, double order);

/* -----------------------
		Matrix class
------------------------*/


class matrix
{
public:
  matrix()
  {
    p = NULL;
    rows = 0;
    cols = 0;
  }

  matrix(const int row_count, const int column_count)
  {
      p = NULL;
      rows = row_count;
      cols = column_count;

      p = new double[rows*cols];
      for (int r = 0; r < rows; r++)
      {

        for (int c = 0; c < cols; c++)
        {
          p[r+c*rows] = 0;
        }
      }
  }

  matrix(const int row_count){
  	rows=row_count;
  	cols=1;

  	p = new double[rows*cols];
  	for (int r = 0; r < rows; r++)
    {
        p[r] = 0;
    }
  }

  matrix(const matrix& a)
  {
    rows = a.rows;
    cols = a.cols;
    p = new double[a.cols*a.rows];

    #pragma omp parallel for
    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        p[r+c*rows] = a.p[r+c*rows];
      }
    }
  }

  matrix(double* q, const int row_count){
    rows = row_count;
    cols = 1;
    p = q;
  }

  matrix(double* q, const int row_count, const int col_count){
    rows = row_count;
    cols = col_count;
    p = q;
  }

  double& operator()(const int r, const int c)
  {
      return p[r+c*rows];
  }

	double& operator()(int i){
				return p[i];
	}

  double get(const int r, const int c) const
  {
      return p[r+c*rows];
  }

  matrix col(const int c)
  {
  	matrix col(rows, 1);

	  #pragma omp parallel for
	  for(int i = 0; i < rows; i++){
	  	col.p[i]=p[i+c*rows];
	  }

    return col;
  }

  matrix subcol(const int c, const int rl, const int rh){
    matrix col(rh-rl, 1);

	  for(int i = 0; i < rh-rl; i++){
	  	col.p[i]=p[i+rl+c*rows];
	  }

    return col;

  }

  void reshape(int new_row, int new_col){
    rows = new_row;
    cols = new_col;
  }

  matrix subvec(const int rl, const int rh){
    matrix col(rh-rl);

	  for(int i = 0; i < rh-rl; i++){
	  	col.p[i]=p[i+rl];
	  }

    return col;
  }

  matrix submatrow(const int rl, const int rh){
    matrix col(rh-rl,cols);
      for(int c = 0; c < cols; c++){
        for(int i = 0; i < rh-rl; i++){
        	col.p[i + c*(rh-rl)]=p[i+rl + c*rows];
        }
    	}
      return col;
  }

  double* point()
  {
  	return p;
  }

  void col_sub(const matrix &col, int c)
  {
     #pragma omp parallel for
  	 for(int i = 0; i < rows; i++){
              p[i+c*rows] = col.p[i];
          }
  }

  matrix row(const int r)
  {
  	matrix row(cols);
	  for(int i = 0; i < cols; i++){
	  	row.p[i]=p[r+rows*i];
	  }
    return row;
  }

  void row_sub(matrix row, int r)
  {
        for(int i = 0; i < cols; i++){
            p[r+i*rows] = row.p[i];
        }
  }

  matrix& operator= (const matrix& a)
  {
    if(rows != a.rows || cols != a.cols || p == NULL){

        if(p != NULL){
            delete [] p;
        }

        rows = a.rows;
        cols = a.cols;

        p = new double[a.cols*a.rows];
    }

    #pragma omp parallel for
    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        p[r+c*rows] = a.p[r+c*rows];
      }
    }
    return *this;
  }

  matrix& Add(const double v)
  {
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r+c*rows] += v;
      }
    }
     return *this;
  }

  matrix& Subtract(const double v)
  {
    return Add(-v);
  }

  matrix& Multiply(const double v)
  {
    #pragma omp parallel for
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r+c*rows] *= v;
      }
    }
     return *this;
  }

  matrix& Divide(const double v)
  {
     return Multiply(1/v);
  }

  friend matrix operator+(const matrix& a, const matrix& b)
  {
      matrix res(a.rows, a.cols);

      #pragma omp parallel for
      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < a.cols; c++)
        {
          res.p[r+c*a.rows] = a.p[r+c*a.rows] + b.p[r+c*a.rows];
        }
      }
      return res;
  }

  friend matrix operator+=(const matrix& a, const matrix& b)
  {
      matrix res(a.rows, a.cols);

      #pragma omp parallel for
      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < a.cols; c++)
        {
          a.p[r+c*a.rows] = a.p[r+c*a.rows] + b.p[r+c*a.rows];
        }
      }
    return matrix();
  }

  friend matrix operator+ (const matrix& a, const double b)
  {
    matrix res = a;
    res.Add(b);
    return res;
  }

  friend matrix operator+ (const double b, const matrix& a)
  {
    matrix res = a;
    res.Add(b);
    return res;
  }


  friend matrix operator- (const matrix& a, const matrix& b)
  {


      matrix res(a.rows, a.cols);

      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < a.cols; c++)
        {
          res.p[r+c*a.rows] = a.p[r+c*a.rows] - b.p[r+c*a.rows];
        }
      }
      return res;
  }


  friend matrix operator- (const matrix& a, const double b)
  {
    matrix res = a;
    res.Subtract(b);
    return res;
  }

  friend matrix operator- (const double b, const matrix& a)
  {
    matrix res = -a;
    res.Add(b);
    return res;
  }

  friend matrix operator- (const matrix& a)
  {
    matrix res(a.rows, a.cols);

    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        res.p[r+c*a.rows] = -a.p[r+c*a.rows];
      }
    }

    return res;
  }


  friend matrix operator* (const matrix& a, const matrix& b)
  {
      matrix res(a.rows, b.cols);
      double sum = 0;

      #pragma omp parallel for private(sum)
      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < b.cols; c++)
        {
          sum = 0;
          for (int c_res = 0; c_res < a.cols; c_res++)
          {
             sum += a.p[r+c_res*a.rows] * b.p[c_res+c*b.rows];
          }
          res.p[r+c*a.rows] = sum;
        }
      }
      return res;
  }

  friend matrix operator* (const matrix& a, const double b)
  {
    matrix res = a;
    res.Multiply(b);
    return res;
  }

  friend matrix operator* (const double b, const matrix& a)
  {
    matrix res = a;
    res.Multiply(b);
    return res;
  }

  friend matrix operator/ (const matrix& a, const matrix& b)
  {
    matrix res(a.rows, a.cols);

    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        res.p[r+c*a.rows] = a.p[r+c*a.rows] / b.p[r+c*a.rows];
      }
    }

    return res;
  }


  friend matrix operator/ (const matrix& a, const double b)
  {
    matrix res = a;
    res.Divide(b);
    return res;
  }


  friend matrix operator/ (const double b, const matrix& a)
  {
    matrix b_matrix(0, 0);
    b_matrix(0,0) = b;

    matrix res = b_matrix / a;
    return res;
  }


  friend matrix operator & (const matrix& a, const matrix& b)
  {
    matrix res(a.rows, a.cols);
    #pragma omp parallel for
    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        res.p[r+c*a.rows] = a.p[r+c*a.rows] * b.p[r+c*a.rows];
      }
    }

    return res;
  }

  matrix Minor(const int row, const int col) const
  {
    matrix res;

    res = matrix(rows - 1, cols - 1);


    for (int r = 0; r < (rows - (row+1 >= rows)); r++)
    {
      for (int c = 0; c < (cols - (col+1 >= cols)); c++)
      {
        res(r - (r > row), c - (c > col)) = p[r+c*rows];
      }
    }

    return res;
  }

  int size(const int i) const
  {
    if (i == 1)
    {
      return rows;
    }
    else if (i == 2)
    {
      return cols;
    }
    return 0;
  }

  int nrows() const
  {
    return rows;
  }

  int ncols() const
  {
    return cols;
  }

  void randu(){
		for(int i = 0; i < rows; i++){
	  		for(int j = 0; j < cols; j++){
				p[i+rows*j]=((double) rand() / (double) (RAND_MAX));
				}
			}
	}

  void zeros(){
		for(int i = 0; i < rows; i++){
	  		for(int j = 0; j < cols; j++){
				p[i+rows*j]=0;
				}
			}
	}

  double determinant(int k)
  {
    double s=1,det=0;
    matrix b(rows, cols);
    int i,j,m,n,c;
    if (k==1)
      {
       return (p[0]);
      }
    else
      {
       det=0;
       for (c=0;c<k;c++)
         {
          m=0;
          n=0;
          for (i=0;i<k;i++)
            {
              for (j=0;j<k;j++)
                {
                  b.p[j*rows+i]=0;
                  if (i != 0 && j != c)
                   {
                     b.p[m*rows+n]=p[j*rows+i];
                     if (n<(k-2))
                      n++;
                     else
                      {
                       n=0;
                       m++;
                       }
                     }
                 }
               }
            det=det + s * (p[c*rows] * b.determinant(k-1));
            s=-1 * s;
            }
      }

      return (det);
  }

  matrix t() const
  {
  	matrix T(cols,rows);
  	for(int i = 0; i < rows; i++){
  		for(int j = 0; j < cols; j++){
  			T(j,i) = p[i+rows*j];
  		}
  	}
    return T;
  }

  matrix cross() const
   {
    double sum;
    matrix res(cols, cols);
    #pragma omp parallel for private(sum)
    for (int r = 0; r < cols; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        sum = 0;
        for (int c_res = 0; c_res < rows; c_res++)
        {
          sum += p[c_res+r*rows] * p[c_res+c*rows];
        }
        res.p[r+c*cols] = sum;
      }
    }


    return res;
  }

  friend matrix operator% (const matrix& a, const matrix& b)
  {
      matrix res(a.cols, b.cols);
      double sum;
      #pragma omp parallel for private(sum)
      for (int r = 0; r < a.cols; r++)
      {
        for (int c = 0; c < b.cols; c++)
        {
          sum = 0;
          for (int c_res = 0; c_res < a.rows; c_res++)
          {
            sum += a.p[c_res+r*a.rows] * b.p[c_res+c*b.rows];
          }
          res.p[r+c*a.cols] = sum;
        }
      }
      return res;
  }

  matrix floorzero() const
  {
    matrix P(rows,cols);
    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        if(p[i+rows*j] > 0){
          P(i,j) = p[i+rows*j];
        }
      }
    }
    return P;
  }

  matrix geq() const
  {
    matrix P(rows,cols);
    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        if(p[i+rows*j] > 0){
          P(i,j) = 1.0;
        }
      }
    }
    return P;
  }

  void print() const
  {
    if (p != NULL)
    {
      printf("[");
      for (int r = 0; r < rows; r++)
      {
        if (r > 0)
        {
          printf(" ");
        }
        for (int c = 0; c < cols-1; c++)
        {
          printf("%.3f, ", p[r+rows*c]);
        }
        if (r < rows-1)
        {
          printf("%.3f;\n", p[r+rows*(cols-1)]);
        }
        else
        {
          printf("%.3f]\n", p[r+rows*(cols-1)]);
        }
      }
    }
    else
    {
      printf("[ ]\n");
    }
  }

  void print(string file){
    ofstream simfile;
    simfile.open(file.c_str());
    simfile.precision(dbl::digits10);
    for(int i = 0; i < rows; i++){
      simfile << p[i];
      for(int j = 1; j < cols; j++){
        simfile << ", " << p[i+rows*j];
      }
      simfile << endl;
    }

    simfile.close();
  }

  void printv() const
  {
    if (p != NULL)
    {
      printf("[");
      for (int r = 0; r < rows; r++)
      {
        if (r > 0)
        {
          printf(" ");
        }
        for (int c = 0; c < cols-1; c++)
        {
          printf("%.3f, ", p[r+rows*c]);
        }
        if (r < rows-1)
        {
          printf("%.3f; ", p[r+rows*(cols-1)]);
        }
        else
        {
          printf("%.3f]", p[r+rows*(cols-1)]);
        }
      }
    }
    else
    {
      printf("[ ]\n");
    }
  }

  matrix cluster(){
     double vi;
     matrix clust(cols, cols);
     for(int a = 0; a < cols; a++){
         for(int b = a; b < cols; b++){
        	    vi = 0;
              for(int i = 0; i < rows; i++){
        	        vi += p[i + a*cols] * p[i + a*cols];
	                for(int j = i + 1; j < rows; j++){
                    vi += 2*p[i + a*cols] * p[j + b*cols];
                  }
              }
              clust(a, b) = vi;
              clust(b, a) = vi;
          }
      }
      return clust;
  }

  virtual ~matrix()
  {
    delete [] p;
    p = NULL;
  }

private:
  int rows;
  int cols;
  double* p;
};

int svd(matrix &mat_arg_a, matrix &mat_arg_w, matrix &mat_arg_v);

int size(const matrix& a, const int i);

matrix ones(const int rows, const int cols);

matrix eye(int rows);

matrix zeros(const int rows, const int cols);
double mean(matrix mat);

double freq(matrix mat);

matrix cmean(matrix mat);

matrix csum(matrix mat);

matrix rsum(matrix mat);
matrix pinv(matrix a);
matrix pinv(double *a, int m, int n);

double var(matrix mat);
double sum(matrix mat);
matrix mabs(matrix mat);
matrix cmeansq(matrix mat);
matrix csqrtmeansq(matrix mat);
matrix sqrt(matrix mat);

matrix cov(matrix mat);
double cov(matrix mat1, matrix mat2);
double cov_vec(matrix mat1, matrix mat2);


matrix diag(const int n);

matrix diag(const matrix& v);

matrix sqrtdiags(const matrix& v);

double det(const matrix& a);

void Swap(double& a, double& b);

matrix inv(matrix a);

matrix chol(matrix A);

matrix sub_f(matrix L, matrix b);

matrix sub_b(matrix U, matrix b);

matrix solve_sym(matrix S, matrix b);

void dsvd(matrix mat, matrix &a, matrix &w, matrix &v);

/* -----------------------
		Integer matrix class
------------------------*/


class imat
{
public:

  imat()
  {


    p = NULL;
    rows = 0;
    cols = 0;
  }


  imat(const int row_count, const int column_count)
  {

    p = NULL;

    if (row_count > 0 && column_count > 0)
    {
      rows = row_count;
      cols = column_count;

      p = new int[rows*cols];
      for (int r = 0; r < rows; r++)
      {

        for (int c = 0; c < cols; c++)
        {
          p[r+c*rows] = 0;
        }
      }
    }
  }

  imat(const int row_count){
  	rows=row_count;
  	cols=1;

  	p = new int[rows*cols];

  	for (int r = 0; r < rows; r++)
    {

        p[r] = 0;
    }
  }


  imat(const imat& a)
  {
    rows = a.rows;
    cols = a.cols;
    p = new int[a.rows*a.cols];
    #pragma omp parallel for
    for (int r = 0; r < a.rows; r++)
    {

      for (int c = 0; c < a.cols; c++)
      {
        p[r+a.rows*c] = a.p[r+a.rows*c];
      }
    }
  }

  int& operator()(const int r, const int c)
  {

      return p[r+rows*c];

  }

	int& operator()(int i){

				return p[i];

	}

  int get(const int r, const int c) const
  {

      return p[r+rows*c];

  }

  imat col(const int c)
  {
  	imat col(rows, 1);

	  for(int i = 0; i < rows; i++){
	  	col.p[i]=p[i+rows*c];
	  }
    return col;

  }

  imat subcol(const int c, const int rl, const int rh){
    imat col(rh-rl, 1);

	  for(int i = 0; i < rh-rl; i++){
	  	col.p[i]=p[i+rl+c*rows];
	  }
    return col;
  }

  imat subvec(const int rl, const int rh){
    imat col(rh-rl);

	  for(int i = 0; i < rh-rl; i++){
	  	col.p[i]=p[i+rl];
	  }

    return col;
  }

  imat row(const int r)
  {
  	imat row(1, cols);

	  for(int i = 0; i < cols; i++){
	  	row.p[i]=p[r+rows*i];
	  }

    return row;

  }


  imat& operator= (const imat& a)
  {
    rows = a.rows;
    cols = a.cols;

    delete [] p;

    p = new int[a.rows*a.cols];
    for (int r = 0; r < a.rows; r++)
    {

      for (int c = 0; c < a.cols; c++)
      {
        p[r+a.rows*c] = a.p[r+a.rows*c];
      }
    }
    return *this;
  }


  imat& Add(const int v)
  {
    #pragma omp parallel for
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r+rows*c] += v;
      }
    }
     return *this;
  }


  imat& Subtract(const int v)
  {
    return Add(-v);
  }


  imat& Multiply(const int v)
  {
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r+rows*c] *= v;
      }
    }
     return *this;
  }


  friend imat operator+ (const imat& a, const int b)
  {
    imat res = a;
    res.Add(b);
    return res;
  }

  friend imat operator+ (const int b, const imat& a)
  {
    imat res = a;
    res.Add(b);
    return res;
  }


  friend imat operator- (const imat& a)
  {
    imat res(a.rows, a.cols);

    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        res.p[r+a.rows*c] = -a.p[r+a.rows*c];
      }
    }

    return res;
  }


  friend imat operator- (const imat& a, const int b)
  {
    imat res = a;
    res.Subtract(b);
    return res;
  }


  friend imat operator- (const int  b, const imat& a)
  {
    imat res = -a;
    res.Add(b);
    return res;
  }

  int size(const int i) const
  {
    if (i == 1)
    {
      return rows;
    }
    else if (i == 2)
    {
      return cols;
    }
    return 0;
  }


  int nrows() const
  {
    return rows;
  }


  int ncols() const
  {
    return cols;
  }

  imat t() const
  {
  	imat T(cols,rows);
  	for(int i = 0; i < cols; i++){
  		for(int j = 0; j < rows; j++){
  			T(j,i) = p[i+rows*j];
  		}
  	}
    return T;
  }


  int* point()
  {
  	return p;
  }


  void print() const
  {
    if (p != NULL)
    {
      printf("[");
      for (int r = 0; r < rows; r++)
      {
        if (r > 0)
        {
          printf(" ");
        }
        for (int c = 0; c < cols-1; c++)
        {
          printf("%i, ", p[r+rows*c]);
        }
        if (r < rows-1)
        {
          printf("%i;\n", p[r+rows*(cols-1)]);
        }
        else
        {
          printf("%i]\n", p[r+rows*(cols-1)]);
        }
      }
    }
    else
    {
      printf("[ ]\n");
    }
  }

  void print(string file){
    ofstream simfile;
    simfile.open(file.c_str());

    for(int i = 0; i < rows; i++){
      simfile << p[i];
      for(int j = 1; j < cols; j++){
        simfile << ", " << p[i+rows*j];
      }
      simfile << endl;
    }

    simfile.close();
  }

public:

  ~imat()
  {

    delete [] p;
    p = NULL;
  }

private:
  int rows;
  int cols;
  int* p;     // pointer to a matrix with doubles
};


matrix exp(matrix mat);

matrix pseudoinv(matrix a);
void svdcmp(matrix &a, matrix &w, matrix &v);
matrix log(matrix mat);

void jacobi_eigenvalue (matrix a, int it_max, matrix &v,
  matrix &d, int &it_num, int &rot_num );

matrix readcsv(string files, int rows, int cols);

struct ord{
	matrix vect;
	imat index;
};

void quickOrderi(double* p1, int* p2, int left, int right);
ord quickOrder(matrix arrr, imat orders, int left, int right);
void quickSort(double* p1, int left, int right);
void quicksort2(double* p1, int left, int right);

#endif
