#include "matrix.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

using namespace std;

#ifndef DAT_H
#define DAT_H
typedef std::numeric_limits< double > dbl;


/*

This is a data class. It's a large matrix with a maximum size where columns can be indexed by name.

Elements
dat       : large matrix which stores the values of all variables with dimensions:
            (nsim - nbr) * nsim X nvar
variables : std::vector of variable names, will be length feed at any given time
nsim      : tracks the total quantity we will simulate
nbr       : tracks burn in quantities
nvar      : (maximum) number of variables
nf        : number of firms/individuals in matrix
feed      : current number of variables taken up
na2       : variable allowing interpolation over profitability shock
zs        : stores random draws to avoid resimulating data (saves time)
is        : stores location of random draws in grid

mul,      : stores previous parameters for AR1 process to see if resimulating
rhol,       is necessary for z draws.
sigmal

zfeed     : checks if all random shocks have been generated
nrep      : number of times data set will be replicated
*/

class data{
	public:
            matrix dat;
            vector<string> variables;
            int nsim, nbr, nvar, feed, na2, nf;
           matrix zs;
            imat   is;
            double mul;
            double rhol;
            double sigmal;
            int zfeed, nrep, resi, seedi;

            data(int nbrl, int nsiml, int nfl, int nvarl, int nrepl = 10, int nstart = 3){
                  nsim = nsiml, nbr = nbrl, nvar = nvarl, nf = nfl;
                  nrep = nrepl;
                  feed = 0;
                  zfeed = 0;
                  resi  = 0;
                  variables.clear();
                  seedi = 145995;
                  dat = matrix((nsim-nbr)*nf, nvar);
                  zs  = matrix(nsim * nf, nrep);
                  is  = imat(nsim * nf, nrep);
            }

            // Clear the data

            data copy(){
              data cdat(nsim, nbr, nvar, nf);
              cdat.feed = feed;
              cdat.variables = variables;
              cdat.zs = zs;
              cdat.is = is;
              cdat.nrep = nrep;
              cdat.zfeed = zfeed;
              cdat.rhol = rhol;
              cdat.seedi = seedi;
              cdat.sigmal = sigmal;
              cdat.mul = mul;
              for(int i = 0; i < nvar; i++){
                for(int j = 0; j < (nsim-nbr)*nf; j++){
                  cdat.dat(j, i) = dat(j, i);
                }
              }

              return cdat;
            }

            void replace_var(string name, matrix varr){
                  int ind = get_varindex(name);
                  dat.col_sub(varr, ind);
            }

            void wipe(){
              feed = 0;
              variables.clear();
            }

            void reset(){
              resi = 0;
	      			zfeed  = 0;
              sigmal = 0;
              rhol   = 0;
              mul    = 0;
              feed = 0;
              variables.clear();
            }

            // Add a variable to the data using a matrix

            void add_variable(string name, matrix &v){

                  variables.push_back(name);

                  if(feed < nvar){
                        dat.col_sub(v, feed);
                        feed++;
                  }
            }
           void add_variable(string name, matrix &vi, matrix &d){
                  matrix v = vi/d;
                  variables.push_back(name);

                  if(feed < nvar){
                        dat.col_sub(v, feed);
                        feed++;
                  }
            }
            // Add a variable to the data using an integer matrix

            void add_variable(string name, imat &v){

                  variables.push_back(name);

                  if(feed < nvar){
                        for(int i = 0; i < nsim - nbr; i++){
                              dat(i, feed) = (double) v(i);
                        }
                        feed++;
                  }
            }

            // Where did we put that variable again?

            int get_varindex(string name){
                  int dex = 100000;
                  for(int i = 0; i < variables.size(); i++){
                        if(variables[i] == name){
                              dex = i;
                        }
                  }



                  return dex;
            }

            int get_varindex_check(string name){
                  int dex = 100000;
                  for(int i = 0; i < variables.size(); i++){
                        if(variables[i] == name){
                              dex = i;
                        }
                  }

                  return dex;
            }

            // Extract the variable as a vector

            matrix get_var(string name){

                  int ind = get_varindex(name);

                  matrix var = dat.col(ind);

                  return var;
            }

            matrix get_var_rm(string name){

                  int ind = get_varindex(name);

                  matrix var = dat.col(ind);
                  matrix varc(var.nrows());

                  double m = mean(var);
                  double msub = 0;

                  for(int i = 0; i < nf; i++){
                    msub = 0;
	                  for(int t = 0; t < nsim - nbr; t++){
                      msub += var(i * (nsim-nbr) + t) / ((double) nsim - nbr);
	                  }

	                  for(int t = 0; t < nsim - nbr; t++){
                      varc(i * (nsim-nbr) + t) = var(i * (nsim-nbr) + t) - msub + m;
	                  }

                  }


                  return varc;
            }

            matrix get_var_dm(string name){

                  int ind = get_varindex(name);

                  matrix var = dat.col(ind);
                  matrix varc(var.nrows());

                  double msub = 0;

                  for(int i = 0; i < nf; i++){
                    msub = 0;
	                  for(int t = 0; t < nsim - nbr; t++){
                      msub += var(i * (nsim-nbr) + t) / ((double) t);
	                  }

	                  for(int t = 0; t < nsim - nbr; t++){
                      varc(i * (nsim-nbr) + t) = var(i * (nsim-nbr) + t) - msub;
	                  }

                  }


                  return varc;
            }

            // Print all this data to a given text file

            void print(string file){
              cout << variables.size() << " " << nvar << " " << dat.nrows() << " " << dat.ncols() << endl;
                ofstream simfile;
                simfile.open(file.c_str());
                simfile.precision(dbl::digits10);
                for(int i = 0; i < variables.size(); i++){
                  simfile << variables[i].c_str() << ", ";
                }



                if(nvar > variables.size()){
                  for(int j = 0; j < nvar - variables.size() - 1; j++){
                    simfile << " ,";
                  }
                }




                simfile << endl;

                for(int i = 0; i < dat.nrows(); i++){
                  simfile << dat(i, 0);
                  for(int j = 1; j < dat.ncols(); j++){
                    simfile << ", " << dat(i, j);
                  }
                  simfile << endl;
                }

                simfile.close();
            }
};

#endif
