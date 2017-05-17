#include "matrix.h"
#include <math.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "epf.h"
#define _USE_MATH_DEFINES

imat permute(int p, int n, int no){
  
  imat options(no, n), retm(p, no);
  for(int i=0; i<no; i++){
    
    for(int j=0; j<n; j++){
      
      options(i,j)=j;
    }
    
    options(i,n-1)=i;
    options(i,i)=n-1;
    
  }
  
	int i_r, t;
	
	for(int i = 0; i < no; i++){
	
		for(int k = 0; k < p; k++){
			
			t = n - 2 - k;
			i_r = ((int) floor((t + 1) * ((double) rand() / (double) RAND_MAX)));
			
			retm(k, i) = options(i, i_r);
			
			options(i, i_r) = options(i, t);
			options(i, t) = retm(k, i);
			
		}
		
	}
	
	return retm;
}

imat randi(int n, int m){
    
	imat ret(m);
	
	for(int i = 0; i < m; i++){
	    
		ret(i)=((int) floor(n * ((double) rand() / (double) RAND_MAX)));
	    
	}
	
	return ret;
}

double cmedian(matrix vecs){
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

int find_nearest(double x, matrix vect){
	int i_u = vect.nrows();
	int i_l = 0;
	int iterl = 0;
	//std::cout << i_u << std::endl;
	if(x <= vect(0)){
		i_l = - 2;
	}
	else if(x >= vect(i_u-1)){
		i_l = vect.nrows() + 2;
	}
	else{
		while(i_u - i_l > 1){
			int i_t = floor((i_u + i_l) / 2);
			//std::cout << i_t << std::endl;
			double xt = vect(i_t);
			if(xt <= x){
				i_l = i_t;
			}
		
			if(xt > x){
				i_u = i_t;
			}
		if(i_u < vect.nrows()){
		    if(vect(i_u) == vect(i_l)){
			    i_u = i_l;
			    printf("UHOH\n");
			    i_l = -2;
		    }
		}
		if(iterl >= vect.nrows()){
			i_u = i_l;
			i_u = - 2;
			i_l = - 2;
		}
		
		iterl++;
		}
		
	}
	return i_l;	
}


void print_mat(matrix mat, string file){
	ofstream simfile;
	simfile.open(file.c_str());
	
	for(int i = 0; i < mat.nrows(); i++){
		simfile << mat(i,0);
		for(int j = 1; j < mat.ncols(); j++){
			simfile << ", " << mat(i,j);
		}
		simfile << endl;
	}
	
	simfile.close();
}

void print_imat(imat mat, string file){
	ofstream simfile;
	simfile.open(file.c_str());
	
	for(int i = 0; i < mat.nrows(); i++){
		simfile << mat(i,0);
		for(int j = 1; j < mat.ncols(); j++){
			simfile << ", " << mat(i,j);
		}
		simfile << endl;
	}
	
	simfile.close();
}


void app_mat(matrix mat, string file, int seed){
	ofstream simfile;
	simfile.open(file.c_str(), ios::app);
	
	simfile << mat(0,0);
	for(int i = 1; i < mat.nrows(); i++){
		simfile << ", " << mat(i,0);
	}
	simfile << ", " << seed << endl;
	
	simfile.close();
}

void app_rmat(matrix mat, string file, int seed, matrix f){
	ofstream simfile;
	simfile.open(file.c_str(), ios::app);
	
	simfile << mat(0,0);
	for(int i = 1; i < mat.nrows(); i++){
		simfile << ", " << mat(i,0);
	}
	
	simfile << ", " << seed;
	
	for(int i = 0; i < f.nrows(); i++){
	    for(int j = 0; j < f.ncols(); j++){
	        simfile << ", " << f(i,j);
	    }
	}
	
	simfile << endl;
	simfile.close();
}
