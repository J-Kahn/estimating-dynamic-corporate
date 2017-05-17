#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
typedef std::numeric_limits< double > dbl;

int size(const matrix& a, const int i)
{
  return a.size(i);
}

int i4_max(int a, int b){
  if(a > b){
    return(a);
  }
  else{
    return(b);
  }
}

int i4_min(int a, int b){
  if(a < b){
    return(a);
  }
  else{
    return(b);
  }
}

double r8_max(double a, double b){
  if(a > b){
    return(a);
  }
  else{
    return(b);
  }
}

double r8_min(double a, double b){
  if(a < b){
    return(a);
  }
  else{
    return(b);
  }
}

double r8_sign(double a){
  if(a < 0){
    return(-1.0);
  }
  else{
    return(1.0);
  }
}

void doubleCopy(double* matrix_old, double* matrix_new, int ns){
  for(int i = 0; i < ns; i++){
    matrix_new[i] = matrix_old[i];
  }
};

void doubleCopy(double* matrix_old, double* matrix_new, int start, int ns){
  for(int i = start; i < ns + start; i++){
    matrix_new[i] = matrix_old[i];
  }
};

void doubleZeros(double* matrix_new, int ns){
  for(int i = 0; i < ns; i++){
    matrix_new[i] = 0;
  }
};

void doubleZeros(double* matrix_new, int start, int ns){
  for(int i = start; i < ns + start; i++){
    matrix_new[i] = 0;
  }
};

void doubleOnes(double* matrix_new, int ns){
  for(int i = 0; i < ns; i++){
    matrix_new[i] = 1.0;
  }
};

void intSeq(int* matrix_new, int ns, int nf){
  for(int i = 0; i < nf -ns; i++){
    matrix_new[i] = ns + i;
  }
};

void doubleAdd(double* matrix_new, double added, int ns){
  for(int i = 0; i < ns; i++){
    matrix_new[i] = matrix_new[i] + added;
  }
};


void doubleAdd(double* matrix_new, double added, int start, int ns){
  for(int i = start; i < ns + start; i++){
    matrix_new[i] = matrix_new[i] + added;
  }
};

void doublePrint(double* p, string file, int rows, int cols){
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
};

void doublePrint(double* p, string file, int rows){
  ofstream simfile;
  simfile.open(file.c_str());
  simfile.precision(dbl::digits10);
  for(int i = 0; i < rows; i++){
    simfile << p[i] << endl;
  }

  simfile.close();
};


void quickOrderi(double* p1, int* p2, int left, int right) {

	int i = left, j = right;
	double tmp;
	int tmp2;
	double pivot = p1[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (p1[i] < pivot)
			i++;
		while (p1[j] > pivot)
			j--;
		if (i <= j) {
			tmp = p1[i];
			tmp2 = p2[i];
			p1[i] = p1[j];
			p1[j] = tmp;
			p2[i] = p2[j];
			p2[j] = tmp2;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j){
		quickOrderi(p1,p2, left, j);
	}
	if (i < right){
		quickOrderi(p1,p2, i, right);
	}


}

ord quickOrder(matrix arrr, imat orders, int left, int right) {
	ord ordr;
	ordr.vect = arrr;
	ordr.index = orders;
	double* p1 = ordr.vect.point();
	int* p2 = ordr.index.point();

	quickOrderi(p1, p2, left, right);


	return ordr;
}

void quickSort(double* p1, int left, int right) {
	int i = left, j = right;
	double tmp;
	double pivot = p1[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (p1[i] < pivot)
			i++;
		while (p1[j] > pivot)
			j--;
		if (i <= j) {
			tmp = p1[i];
			p1[i] = p1[j];
			p1[j] = tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j){
		quickSort(p1, left, j);
	}
	if (i < right){
		quickSort(p1, i, right);
	}
}


/**
 * returns a matrix with size cols x rows with ones as values
 */
matrix ones(const int rows, const int cols)
{
  matrix res = matrix(rows, cols);

  for (int r = 0; r < rows; r++)
  {
    for (int c = 0; c < cols; c++)
    {
      res(r, c) = 1;
    }
  }
  return res;
}

/*
  Takes absolute value of each item in matrix.
*/
matrix mabs(matrix mat)
{
  matrix res = matrix(mat.nrows(), mat.ncols());

  for (int r = 0; r < mat.nrows(); r++)
  {
    for (int c = 0; c < mat.ncols(); c++)
    {
      res(r, c) = std::abs(mat(r,c));
    }
  }
  return res;
}

/**
 * returns a matrix with size cols x rows with zeros as values
 */
matrix zeros(const int rows, const int cols)
{
  return matrix(rows, cols);
}

double mean(matrix mat){
	double m = 0;

	#pragma omp parallel for reduction(+:m)
	for(int i = 0; i < mat.nrows(); i++){
		m+=mat(i);
	}

	m/=((double) mat.nrows());

	return m;
}

double central_moment(matrix mat, double order){
  double m = mean(mat);
  double mom = 0;

  for(int i = 0; i < mat.nrows(); i++){
    mom += std::pow(mat(i) - m, order);
  }

  mom/=((double) mat.nrows());

  return mom;
}

matrix central_moment_infl(matrix mat, double order){
  double m = mean(mat);
  double c = central_moment(mat, order);
  matrix infl(mat.nrows());

  for(int i = 0; i < mat.nrows(); i++){
    infl(i) = std::pow(mat(i) - m, order) - c;
  }

  return infl;
}

double hp_autocorr(matrix mat, matrix matl){
  double denom = 0;
  double num   = 0;
  double ml;

  for(int i = 0; i < mat.nrows(); i++){
    ml = matl(i);
    denom += ml * mat(i);
    num   += ml * ml;
  }

  return 2 * denom / num + 1.0;
}

matrix hp_autocorr_infl(matrix mat, matrix matl){
  double denom = 0;
  double num   = 0;
  double ml;

  int n = mat.nrows();

  matrix infl(n);

  for(int i = 0; i < n; i++){
    ml      = matl(i);
    denom += ml * mat(i);
    num   += ml * ml;
  }

  double beta =  2 * denom / num + 1.0;

  for(int i = 0; i < n; i++){
    ml      = matl(i);
    infl(i) = - ml / num * (2 * mat(i) + ml - beta * ml) * ((double) n);
  }

  return infl;
}



double freq(matrix mat){
	double m = 0;

	for (int i = 0; i < mat.nrows(); i++){
		if (mat(i) > 0){
			m += 1.0;
		}
	}

	m /= ((double)mat.nrows());

	return m;
}


double sum(matrix mat){
	double m = 0;

	for(int i = 0; i < mat.nrows(); i++){
		m+=mat(i);
	}

	return m;
}

double median(matrix vecs){
  matrix vect=vecs;
  int nt = vect.nrows();
  quickSort(vect.point(),0,nt-1);
  double med;

  if(nt%2==0){
    med=0.5*(vect(nt/2-1)+vect(nt/2));
  }
  else{
    med=vect((nt-1)/2);
  }

  return med;
}


matrix cmean(matrix mat){
	matrix m(mat.ncols());
	#pragma omp parallel for
	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			m(i)+=mat(j,i);
		}
	}

	return m/((double) mat.nrows());
}

matrix cmeansq(matrix mat){
	matrix m(mat.ncols());

	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			m(i)+=mat(j,i)*mat(j,i);
		}
	}

	return m/((double) mat.nrows());
}

matrix sqrt(matrix mat){
  matrix ret(mat.nrows(), mat.ncols());
  #pragma omp parallel for
	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			ret(j, i) = sqrt(mat(j,i));
		}
	}

  return ret;
}

matrix csqrtmeansq(matrix mat){
	matrix m(mat.ncols());

	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			m(i)+=mat(j,i)*mat(j,i);
		}
		m(i) = sqrt(m(i));
	}

	return m/((double) mat.nrows());
}



matrix csum(matrix mat){
	matrix m(mat.ncols());

	#pragma omp parallel for
	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			m(i)+=mat(j,i);
		}
	}

	return m;
}

matrix rsum(matrix mat){
	matrix m(mat.nrows());

	for(int i = 0; i < mat.ncols(); i++){
		for(int j = 0; j < mat.nrows(); j++){
			m(j)+=mat(j,i);
		}
	}

	return m;
}

double var(matrix mat){
	double m = mean(mat);
	double v = 0;
        double z = 0;
	#pragma omp parallel for private(z) reduction(+:v)
	for(int i = 0; i < mat.nrows(); i++){
		z = mat(i) - m;
		v += z * z;
	}

	v /= ((double) mat.nrows());

	return v;
}

double cov_vec(matrix mat1, matrix mat2){
	double m1 = mean(mat1);
	double m2 = mean(mat2);
	double v = 0;
        double z = 0;
	#pragma omp parallel for private(z) reduction(+:v)
	for(int i = 0; i < mat1.nrows(); i++){
		v += (mat1(i) - m1) * (mat2(i) - m2);
	}

	v /= ((double) mat1.nrows());

	return v;
}

matrix cov(matrix mat){

	matrix dm = zeros(mat.nrows(), mat.ncols());
	matrix m = cmean(mat);

	#pragma omp parallel for
	for(int i = 0; i < mat.nrows(); i++){
		for(int j = 0; j < mat.ncols(); j++){
			dm(i, j) = mat(i, j) - m(j);
		}
	}

	matrix cv = (dm.t()*dm)/((double) dm.nrows());

	return cv;
}


double cov(matrix mat1, matrix mat2){

  matrix dm1(mat1.nrows()), dm2(mat2.nrows());
  double m1 = mean(mat1), m2 = mean(mat2);

  #pragma omp parallel for
  for(int i = 0; i < mat1.nrows(); i++){
    dm1(i) = mat1(i) - m1;
    dm2(i) = mat2(i) - m2;
  }

  matrix cv = (dm1.cross())/((double) dm1.nrows());
  double ccv = cv(0);

  return ccv;
}


matrix diag(const int n)
{
  matrix res(n, n);
  for (int i = 0; i < n; i++)
  {
    res(i, i) = 1;
  }
  return res;
}


matrix diag(const matrix& v)
{
  matrix res;
  if (v.ncols() == 1){

    int rows = v.nrows();
    res = matrix(rows, rows);


    for (int r=0; r < rows; r++)
    {
      res(r, r) = v.get(r, 0);
    }
  }
  else if (v.nrows() == 0)
  {

    int cols = v.ncols();
    res = matrix(cols, cols);


    for (int c=0; c < cols; c++)
    {
      res(c, c) = v.get(0, c);
    }
  }

  return res;
}

matrix exp(matrix mat)
{
    matrix expm(mat.nrows(), mat.ncols());
    for(int i = 0; i < mat.nrows(); i++){
        for(int j = 0; j < mat.ncols(); j++){
            expm(i, j) = exp(mat(i, j));
        }
    }

    return expm;
}

matrix log(matrix mat)
{
    matrix expm(mat.nrows(), mat.ncols());
    for(int i = 0; i < mat.nrows(); i++){
        for(int j = 0; j < mat.ncols(); j++){
            expm(i, j) = log(mat(i, j));
        }
    }

    return expm;
}


matrix diags(const matrix& v)
{
  matrix res(v.nrows());

  for (int c=0; c < v.nrows(); c++)
  {
    res(c) = v.get(c,c);
  }

  return res;
}

matrix sqrtdiags(const matrix& v)
{
  matrix res(v.nrows());
  for (int c=0; c < v.nrows(); c++)
  {
    res(c) = sqrt(v.get(c,c));
  }

  return res;
}


/*
 * returns the determinant of matrix a
 */
double det(const matrix& a)
{
  double d = 0;    // value of the determinant
  int rows = a.nrows();
  int cols = a.ncols();

    if (rows == 1)
    {
      d = a.get(0, 0);
    }
    else if (rows == 2)
    {
      d = a.get(0, 0) * a.get(1, 1) - a.get(1, 0) * a.get(0, 1);
    }
    else
    {
      for (int c = 0; c < cols; c++)
      {
        matrix M = a.Minor(0, c);

        d += ((c+1)%2 + (c+1)%2 - 1) * a.get(0, c) * det(M); // faster than with pow()

      }
    }
  return d;
}

matrix eye(int rows){
    matrix I = zeros(rows, rows);
    for(int i = 0; i < rows; i++){
        I(i,i) = 1;
    }
    return I;
}

matrix inv(matrix B){
    matrix A = B;
    int rows = A.nrows();
    matrix I = eye(rows);
    double ratio;

    for(int z = 0; z < rows; z++){
        A = B;
        for (int i = 0; i < rows; i++) {
           for (int j = i+1; j < rows; j++) {
               ratio = A(j,i)/A(i,i) ;
               for (int k = i; k < rows; k++) {
                    A(j,k) -= (ratio*A(i,k));
               }
               I(j,z) -= (ratio*I(i,z));
            }

            I(i, z) = I(i, z)/A(i,i);
            for(int k = i+1; k < rows; k++){
                A(i, k) = A(i,k)/A(i,i);
            }

            A(i, i) = 1;
        }

        for(int i = rows - 1; i >= 0; i--){
            for(int j = i - 1; j >= 0; j--){
               ratio = A(j,i);
               for (int k = i; k < rows; k++) {
                    A(j,k) -= (ratio*A(i,k));
               }
               I(j,z) -= (ratio*I(i,z));
            }
        }
    }

    return I;
}


void swap(int i,int j, double *a){
    int temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}


void quicksort2(double *arr, int left, int right){
    int min = (left+right)/2;
    cout<<"QS:"<<left<<","<<right<<"\n";

    int i = left;
    int j = right;
    double pivot = arr[min];

    while(left<j || i<right)
    {
        while(arr[i]<pivot)
        i++;
        while(arr[j]>pivot)
        j--;

        if(i<=j){
            swap(i,j,arr);
            i++;
            j--;
        }
        else{
            if(left<j)
                quicksort2(arr, left, j);
            if(i<right)
                quicksort2(arr,i,right);
            return;
        }
    }
}

matrix chol(matrix A){

	int n = A.nrows();

	matrix L(n,n);

	L=zeros(n,n);

	for(int j = 0; j < n; j++){

		for(int i = 0; i < j + 1; i++)
		{

			for(int k = 0; k < i; k++)
			{

				L(j, i) += L(i,k) * L(j,k);

			}

			if(i == j)
			{

				L(j, i) = sqrt(A(j, i) - L(j, i));

			}
			else{
				L(j, i) = 1.0/ L(i, i) * (A(j, i) - L(j, i));
			}
		}
	}

	return L;
}

matrix sub_f(matrix L, matrix b){
	int n =  b.nrows(), k = b.ncols();

	matrix x=zeros(n,k);
	for(int l = 0; l < k; l++){
  	for(int i = 0; i < n; i++)
  	{
  		x(i,l) = b(i,l)/L(i, i);

  		for(int j = i; j < n; j++)
  		{
  			b(j, l) = b(j, l) - L(j, i) * x(i, l);
  		}
  	}
  }

	return x;
}

matrix sub_b(matrix U, matrix b){
  int n =  b.nrows(), k = b.ncols();

	matrix x=zeros(n,k);

  for(int l = 0; l < k; l++){
  	for(int i = n - 1; i >= 0; i--)
  	{
  		x(i, l) = b(i, l)/U(i, i);

  		for(int j = 0; j < i; j++)
  		{
  			b(j, l) = b(j, l) - U(j, i) * x(i, l);
  		}
  	}
  }

	return x;
}

matrix solve_sym(matrix S, matrix b){
	matrix L = chol(S);
	matrix y = sub_f(L, b);
	matrix x = sub_b(L.t(),y);
	return x;
}

matrix solve_sym_mult(matrix S, matrix b){
	matrix x(S.nrows(), b.ncols());
	for(int i = 0; i < b.ncols(); i++){
	    x.col_sub(solve_sym(S, b.col(i)), i);
	}
        return x;
}

matrix readcsv(string files, int rows, int cols){
    ifstream file(files.c_str());
    matrix datas(rows,cols);
    double d1;
    for(int row = 0; row < rows; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);
        for (int col = 0; col < cols; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            std::stringstream convertor(val);
            convertor >> datas(row,col);
        }
    }
    return datas;
}






void dscal ( int n, double sa, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}




double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}












double dnrm2 ( int n, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
  double absxi;
  int i;
  int ix;
  double norm;
  double scale;
  double ssq;
  double value;
  cout << n << " " << incx << " ... ";
  if ( n < 1 || incx < 1 )
  {
    norm = 0.0;
    //cout << 0.0 << endl;
  }
  else if ( n == 1 )
  {
    norm = fabs ( x[0] );
   // cout << x[0] << endl;
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
   //   cout << x[i] << " - ";
      if ( x[ix] != 0.0 )
      {
        absxi = fabs ( x[ix] );
        if ( scale < absxi )
        {
          ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
          scale = absxi;
        }
        else
        {
          ssq = ssq + ( absxi / scale ) * ( absxi / scale );
        }
      }
 //     cout << ssq << ". ";
      ix = ix + incx;
    }

    norm  = scale * sqrt ( ssq );
  }
 // cout << endl;
  return norm;
}






void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}






void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s )

//****************************************************************************80
//
//  Purpose:
//
//    DROT applies a plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, double C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
  int i;
  int ix;
  int iy;
  double stemp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[i] + s * y[i];
      y[i]  = c * y[i] - s * x[i];
      x[i]  = stemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[ix] + s * y[iy];
      y[iy] = c * y[iy] - s * x[ix];
      x[ix] = stemp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}






void dswap ( int n, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;
  double temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}





void drotg ( double *sa, double *sb, double *c, double *s )

//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
  double r;
  double roe;
  double scale;
  double z;

  if ( fabs ( *sb ) < fabs ( *sa ) )
  {
    roe = *sa;
  }
  else
  {
    roe = *sb;
  }

  scale = fabs ( *sa ) + fabs ( *sb );

  if ( scale == 0.0 )
  {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
  }
  else
  {
    r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
                     + ( *sb / scale ) * ( *sb / scale ) );
    r = r8_sign ( roe ) * r;
    *c = *sa / r;
    *s = *sb / r;
  }

  if ( 0.0 < fabs ( *c ) && fabs ( *c ) <= *s )
  {
    z = 1.0 / *c;
  }
  else
  {
    z = *s;
  }

  *sa = r;
  *sb = z;

  return;
}






int dsvdc ( double a[], int lda, int m, int n, double s[], double e[],
  double u[], int ldu, double v[], int ldv, double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{

  /*cout << "INPUTS " << endl;
  cout << lda << " " << m << " " << n << " " << ldu << " " << ldv << endl << endl;
  cout << "MATRICES" << endl;
  for(int i = 0; i < n; i ++){
    cout << work[i] << ": ";
    for(int j = 0; j < m; j++){
      cout << a[i + j * m] << " " << (a + j*m)[i] << endl;
    }
    cout << endl;
  }*/

  double b;
  double c;
  double cs;
  double el;
  double emm1;
  double f;
  double g;
  int i;
  int info;
  int iter;
  int j;
  int jobu;
  int k;
  int kase;
  int kk;
  int l;
  int ll;
  int lls;
  int ls;
  int lu;
  int maxit = 30;
  int mm;
  int mm1;
  int mn;
  int mp1;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  double t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  double ztest;
//
//  Determine what is to be computed.
//
  /*for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      cout << a[i + j * m];
    }
    cout << endl;
  }*/
  info = 0;
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
  {
    ncu = i4_min ( m, n );
  }
  else
  {
    ncu = m;
  }

  if ( jobu != 0 )
  {
    wantu = true;
  }

  if ( ( job % 10 ) != 0 )
  {
    wantv = true;
  }
//
//  Reduce A to bidiagonal form, storing the diagonal elements
//  in S and the super-diagonal elements in E.
//
  nct = i4_min ( m-1, n );
  nrt = i4_max ( 0, i4_min ( m, n-2 ) );
  lu = i4_max ( nct, nrt );
  //cout << "nct " << nct << " nrt " << nrt << " lu " << lu << endl;

  for ( l = 1; l <= lu; l++ )
  {
//
//  Compute the transformation for the L-th column and
//  place the L-th diagonal in S(L).
//
    if ( l <= nct )
    {
      s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );
      //cout << " s " << s[l-1] << " " << a[0] <<   endl;
      if ( s[l-1] != 0.0 )
      {
        if ( a[l-1+(l-1)*lda] != 0.0 )
        {
          s[l-1] = r8_sign ( a[l-1+(l-1)*lda] ) * fabs ( s[l-1] );
        }
        dscal ( m-l+1, 1.0 / s[l-1], a+l-1+(l-1)*lda, 1 );
        a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
      }
      s[l-1] = -s[l-1];
      //cout << " sagain " << s[l-1] << endl;
    }

    for ( j = l+1; j <= n; j++ )
    {
//
//  Apply the transformation.
//
      if ( l <= nct && s[l-1] != 0.0 )
      {
        t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 )
          / a[l-1+(l-1)*lda];
        daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
      }
    //cout << "t " << t << endl;
//  Place the L-th row of A into E for the
//  subsequent calculation of the row transformation.
//
      e[j-1] = a[l-1+(j-1)*lda];
      //cout << "ej " << e[j-1];
    }
    //cout << endl;
//
//  Place the transformation in U for subsequent back multiplication.
//
    if ( wantu && l <= nct )
    {
      for ( i = l; i <= m; i++ )
      {
        u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
      }
    }

    if ( l <= nrt )
    {
//
//  Compute the L-th row transformation and place the
//  L-th superdiagonal in E(L).
//
      e[l-1] = dnrm2 ( n-l, e+l, 1 );

      if ( e[l-1] != 0.0 )
      {
        if ( e[l] != 0.0 )
        {
          e[l-1] = r8_sign ( e[l] ) * fabs ( e[l-1] );
        }
        dscal ( n-l, 1.0 / e[l-1], e+l, 1 );
        e[l] = 1.0 + e[l];
      }

      e[l-1] = -e[l-1];
//
//  Apply the transformation.
//
      if ( l+1 <= m && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= m; j++ )
        {
          work[j-1] = 0.0;
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );
        }
       // cout << "workj " << work[j - 1] << endl;
      }
//
//  Place the transformation in V for subsequent back multiplication.
//
      if ( wantv )
      {
        for ( j = l+1; j <= n; j++ )
        {
          v[j-1+(l-1)*ldv] = e[j-1];
        }
      }
    }
  }
//
//  Set up the final bidiagonal matrix of order MN.
//
  mn = i4_min ( m + 1, n );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < n )
  {
    s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
  }

  if ( m < mn )
  {
    s[mn-1] = 0.0;
  }

  if ( nrtp1 < mn )
  {
    e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
  }

  e[mn-1] = 0.0;
//
//  If required, generate U.
//
  if ( wantu )
  {
    for ( i = 1; i <= m; i++ )
    {
      for ( j = nctp1; j <= ncu; j++ )
      {
        u[(i-1)+(j-1)*ldu] = 0.0;
      }
    }

    for ( j = nctp1; j <= ncu; j++ )
    {
      u[j-1+(j-1)*ldu] = 1.0;
    }

    for ( ll = 1; ll <= nct; ll++ )
    {
      l = nct - ll + 1;

      if ( s[l-1] != 0.0 )
      {
        for ( j = l+1; j <= ncu; j++ )
        {
          t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 )
            / u[l-1+(l-1)*ldu];
          daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
        }

        dscal ( m-l+1, -1.0, u+(l-1)+(l-1)*ldu, 1 );
        u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
        for ( i = 1; i <= l-1; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        //cout << "uthing " << u[l-1+(l-1)*ldu];

        }
      }
      else
      {
        for ( i = 1; i <= m; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
        u[l-1+(l-1)*ldu] = 1.0;
      //  cout << "uthing " << u[l-1+(l-1)*ldu];
      }
      //cout << endl;

    }
  }
//
//  If it is required, generate V.
//
  if ( wantv )
  {
    for ( ll = 1; ll <= n; ll++ )
    {
      l = n - ll + 1;

      if ( l <= nrt && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= n; j++ )
        {
          t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 )
            / v[l+(l-1)*ldv];
          daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
        }

      }
      for ( i = 1; i <= n; i++ )
      {
        v[i-1+(l-1)*ldv] = 0.0;
      }
      v[l-1+(l-1)*ldv] = 1.0;
    }
  }
//
//  Main iteration loop for the singular values.
//
  mm = mn;
  iter = 0;

  while ( 0 < mn )
  {
//
//  If too many iterations have been performed, set flag and return.
//
    if ( maxit <= iter )
    {
      info = mn;
      return info;
    }
//
//  This section of the program inspects for
//  negligible elements in the S and E arrays.
//
//  On completion the variables KASE and L are set as follows:
//
//  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
//  KASE = 2     if S(L) is negligible and L < MN
//  KASE = 3     if E(L-1) is negligible, L < MN, and
//               S(L), ..., S(MN) are not negligible (QR step).
//  KASE = 4     if E(MN-1) is negligible (convergence).
//
    for ( ll = 1; ll <= mn; ll++ )
    {
      l = mn - ll;

      if ( l == 0 )
      {
        break;
      }

      test = fabs ( s[l-1] ) + fabs ( s[l] );
      ztest = test + fabs ( e[l-1] );

      if ( ztest == test )
      {
        e[l-1] = 0.0;
        break;
      }
    }

    if ( l == mn - 1 )
    {
      kase = 4;
    }
    else
    {
      mp1 = mn + 1;

      for ( lls = l+1; lls <= mn+1; lls++ )
      {
        ls = mn - lls + l + 1;

        if ( ls == l )
        {
          break;
        }

        test = 0.0;
        if ( ls != mn )
        {
          test = test + fabs ( e[ls-1] );
        }

        if ( ls != l + 1 )
        {
          test = test + fabs ( e[ls-2] );
        }

        ztest = test + fabs ( s[ls-1] );

        if ( ztest == test )
        {
          s[ls-1] = 0.0;
          break;
        }

      }

      if ( ls == l )
      {
        kase = 3;
      }
      else if ( ls == mn )
      {
        kase = 1;
      }
      else
      {
        kase = 2;
        l = ls;
      }
    }

    l = l + 1;
//
//  Deflate negligible S(MN).
//
    if ( kase == 1 )
    {
      mm1 = mn - 1;
      f = e[mn-2];
      e[mn-2] = 0.0;

      for ( kk = 1; kk <= mm1; kk++ )
      {
        k = mm1 - kk + l;
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;

        if ( k != l )
        {
          f = -sn * e[k-2];
          e[k-2] = cs * e[k-2];
        }

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );
        }
      }
    }
//
//  Split at negligible S(L).
//
    else if ( kase == 2 )
    {
      f = e[l-2];
      e[l-2] = 0.0;

      for ( k = l; k <= mn; k++ )
      {
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;
        f = - sn * e[k-1];
        e[k-1] = cs * e[k-1];
        if ( wantu )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
        }
      }
    }
//
//  Perform one QR step.
//
    else if ( kase == 3 )
    {
//
//  Calculate the shift.
//
      scale = r8_max ( fabs ( s[mn-1] ),
              r8_max ( fabs ( s[mn-2] ),
              r8_max ( fabs ( e[mn-2] ),
              r8_max ( fabs ( s[l-1] ), fabs ( e[l-1] ) ) ) ) );

      sm = s[mn-1] / scale;
      smm1 = s[mn-2] / scale;
      emm1 = e[mn-2] / scale;
      sl = s[l-1] / scale;
      el = e[l-1] / scale;
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
      c = ( sm * emm1 ) * ( sm * emm1 );
      shift = 0.0;

      if ( b != 0.0 || c != 0.0 )
      {
        shift = sqrt ( b * b + c );
        if ( b < 0.0 )
        {
          shift = -shift;
        }
        shift = c / ( b + shift );
      }

      f = ( sl + sm ) * ( sl - sm ) - shift;
      g = sl * el;
//
//  Chase zeros.
//
      mm1 = mn - 1;

      for ( k = l; k <= mm1; k++ )
      {
        drotg ( &f, &g, &cs, &sn );

        if ( k != l )
        {
          e[k-2] = f;
        }

        f = cs * s[k-1] + sn * e[k-1];
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * s[k];
        s[k] = cs * s[k];

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
        }

        drotg ( &f, &g, &cs, &sn );
        s[k-1] = f;
        f = cs * e[k-1] + sn * s[k];
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];

        if ( wantu && k < m )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
        }
      }
      e[mn-2] = f;
      iter = iter + 1;
    }
//
//  Convergence.
//
    else if ( kase == 4 )
    {
//
//  Make the singular value nonnegative.
//
      if ( s[l-1] < 0.0 )
      {
        s[l-1] = -s[l-1];
        if ( wantv )
        {
          dscal ( n, -1.0, v+0+(l-1)*ldv, 1 );
        }
      }
//
//  Order the singular value.
//
      for ( ; ; )
      {
        if ( l == mm )
        {
          break;
        }

        if ( s[l] <= s[l-1] )
        {
          break;
        }

        t = s[l-1];
        s[l-1] = s[l];
        s[l] = t;

        if ( wantv && l < n )
        {
          dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
        }

        if ( wantu && l < m )
        {
          dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
        }

        l = l + 1;
      }
      iter = 0;
      mn = mn - 1;
    }
  }

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
     // cout << u[i + j * n] << " ";
    }
    //cout << endl;
  }

  return info;
}



void r8mat_svd_linpack ( int m, int n, double a[], double u[], double s[],
  double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
//
//  Discussion:
//
//    The singular value decomposition of a real MxN matrix A has the form:
//
//      A = U * S * V'
//
//    where
//
//      U is MxM orthogonal,
//      S is MxN, and entirely zero except for the diagonal;
//      V is NxN orthogonal.
//
//    Moreover, the nonzero entries of S are positive, and appear
//    in order, from largest magnitude to smallest.
//
//    This routine calls the LINPACK routine DSVDC to compute the
//    factorization.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix A.
//
//    Input, double A[M*N], the matrix whose singular value
//    decomposition we are investigating.
//
//    Output, double U[M*M], S[M*N], V[N*N], the factors
//    that form the singular value decomposition of A.
//
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int lwork;
  double *sdiag;
  double *work;
//
//  The correct size of E and SDIAG is min ( m+1, n).
//
  a_copy = new double[m*n];
  e = new double[m*n];
  sdiag = new double[m*n];
  work = new double[m*m];
//
//  Compute the eigenvalues and eigenvectors.
//
  job = 11;
  lda = m;
  ldu = m;
  ldv = n;
//
//  The input matrix is destroyed by the routine.  Since we need to keep
//  it around, we only pass a copy to the routine.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a_copy[i+j*m] = a[i+j*m];
     // cout << a_copy[i + j *m];
    }
   // cout << endl;
  }
  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "R8MAT_SVD_LINPACK - Failure!\n";
    cout << "  The SVD could not be calculated.\n";
    cout << "  LINPACK routine DSVDC returned a nonzero\n";
    cout << "  value of the error flag, INFO = " << info << "\n";
    return;
  }
//
//  Make the MxN matrix S from the diagonal values in SDIAG.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        s[i+j*m] = sdiag[i];
      }
      else
      {
        s[i+j*m] = 0.0;
      }
    }
  }
//
//  Note that we do NOT need to transpose the V that comes out of LINPACK!
//
 /* delete [] a_copy;
  delete [] e;
  delete [] sdiag;
  delete [] work;*/

  return;
}


double *pseudo_inverse ( int m, int n, double u[], double s[],
  double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PSEUDO_INVERSE computes the pseudoinverse.
//
//  Discussion:
//
//    Given the singular value decomposition of a real MxN matrix A:
//
//      A = U * S * V'
//
//    where
//
//      U is MxM orthogonal,
//      S is MxN, and entirely zero except for the diagonal;
//      V is NxN orthogonal.
//
//    the pseudo inverse is the NxM matrix A+ with the form
//
//      A+ = V * S+ * U'
//
//    where
//
//      S+ is the NxM matrix whose nonzero diagonal elements are
//      the inverses of the corresponding diagonal elements of S.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix A.
//
//    Input, double A[M*N], the matrix whose singular value
//    decomposition we are investigating.
//
//    Input, double U[M*M], S[M*N], V[N*N], the factors
//    that form the singular value decomposition of A.
//
//    Output, double PSEUDO_INVERSE[N*M], the pseudo_inverse of A.
//
{
  double *a_pseudo;
  int i;
  int j;
  int k;
  double *sp;
  double *sput;

  sp = new double[n*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j && s[i+i*m] != 0 )
      {
        sp[i+j*n] = 1.0 / s[i+i*m];
      }
      else
      {
        sp[i+j*n] = 0.0;
      }
    }
  }

  sput = new double[n*m];
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      sput[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        sput[i+j*n] = sput[i+j*n] + sp[i+k*n] * u[j+k*m];
      }
    }
  }

  a_pseudo = new double[n*m];
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      a_pseudo[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        a_pseudo[i+j*n] = a_pseudo[i+j*n] + v[i+k*n] * sput[k+j*n];
      }
    }
  }

  //delete [] sp;

  return a_pseudo;
}


matrix pinv(matrix a)
{
  int n = a.ncols();
  int m = a.nrows();
  double *a_pseudo_d;
  matrix s(m,n), u(m,m), v(n,n);

  /*for(int i = 0; i < m; i ++){
    for(int j = 0; j < m; j++){
      cout << a.point()[j + i * m] << " ";
    }
    cout << endl;
  }*/

  r8mat_svd_linpack( m, n, a.point(), u.point(), s.point(), v.point() );
  u.print();
  v.print();
  s.print();
  //v = v.t();

  a_pseudo_d = pseudo_inverse(m, n, u.point(), s.point(), v.point());

  matrix a_pseudo(n,m);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      a_pseudo(i,j) = a_pseudo_d[i + j * n];
    }
  }

 // (a * a_pseudo).print();
  return a_pseudo;
}

matrix pinv(double *a, int n, int m)
{
  double *a_pseudo_d;
  matrix s(m,n), u(m,m), v(n,n);

  for(int i = 0; i < m; i ++){
    for(int j = 0; j < m; j++){
      cout << a[j + i * m] << " ";
    }
    cout << endl;
  }

  r8mat_svd_linpack ( m, n, a, u.point(), s.point(), v.point() );
  u.print();
  v.print();
  s.print();
  //v = v.t();

  a_pseudo_d = pseudo_inverse(m, n, u.point(), s.point(), v.point());

  matrix a_pseudo(n,m);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      a_pseudo(i,j) = a_pseudo_d[i + j * n];
    }
  }

  //(a * a_pseudo).print();
  return a_pseudo;
}
