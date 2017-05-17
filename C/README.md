

# Usage

 Here's how to make an executable for a normal monte carlo or simulation:

	make monte monte=mc_epfq_drs.cpp optim=simann.cpp simulate=simulate_drs.cpp model=model_drs.cpp

This compiles to an executable "monte" using optimization from simann.cpp and running the estimation in mc_epfq_drs.cpp, the drs model in model_drs.cpp and the simulation routine in simulate_drs_quick.cpp. This mostly serves to reduce the number of files hanging around in this folder. The same procedure works for compiling files which simulate the data for the Monte Carlo trials, for instance:

	make monte monte=mc_sim_drs.cpp optim=simann.cpp simulate=simulate_drs.cpp model=model_drs.cpp

compiles an executable "monte" that will simulate the drs model.

# Files

+ regress.cpp  : This file has all the code behind running linear regressions and setting up polynomial EPFs.

+ data.h           : This has a class which is used for storing data from simulations with general

+ epf.cpp          : Code for calculating empirical policy functions and moments for simulated data.

+ epf.h            : Header file with references to epf.cpp.

+ inference.cpp    :  Code to calculate jacobians for inference

+ matrix.h         :  Header file for matrices.

+ matrixu.cpp : Utility file for matrices.

+ neldmead.cpp     : Contains code for Nelder-Mead optimization, currently unused.

+ objec.cpp        : Objective files for various types of estimation.

+ opt.h            : Header file for simann, neldmead and objec functions.

+ simann.cpp       : Contains code for simulated annealing optimization.

+ simulate_{MODEL}.cpp     : Contains the code to simulate DRS or CRS model.

+ simulation.h     : Header file to go with simulation.cpp.

+ solve.h      : Contains class dynprob for setting up and solving dynamic problems.

+ utilities.cpp    : A bunch of random functions.

+ mc_{   }.cpp : Runs the monte carlo for a particular model or estimator.

+ mc_sim_{    }.cpp : Simulates the data for Monte Carlos for a particular estimator.

# Objects

+ dynprob : Class contains policy functions, flow payoffs and vfi() function for solving dynamic problems. See solve.h.

+ data    : Data class has a large matrix for storing values from a simulation, indexed by name so I can easily add variables to be simulated. See data.h.

+ matrix  : Standard matrix class. See matrix.h.

+ imat    : Integer matrix class. See matrix.h.

# Functions

## OPTIMIZATION

+ simann   : Simulated annealing code. Finds minimum of objective function, objec, of the form in objec.cpp. Ported from Bill Goffe's Fortran code.

	Arguments:

	- N               : Number of variables in the function to be optimized. (INT)

	- X               : The starting values for the variables of the function to be optimized. (DP(N))

	- MAX             : Denotes whether the function should be maximized or minimized. A true value denotes maximization while a false value denotes minimization. (L)

	- RT              : The temperature reduction factor.  The value suggested by Corana et al. is .85. See Goffe et al. for more advice. (DP)

	- EPS             : Error tolerance for termination. (EP)

	- NS              : Number of cycles.  After NS*N function evaluations, each element of VM is adjusted so that approximately half of all function evaluations are accepted.  The suggested value is 20. (INT)

	- NT              : Number of iterations before temperature reduction.

	- NEPS            : Number of final function values used to decide upon termination.  See EPS.  Suggested value is 4. (INT)

	- MAXEVL          : The maximum number of function evaluations.  If it is exceeded, IER = 1. (INT)

	- LB              : The lower bound for the allowable solution variables. (DP(N))

	- UB              : The upper bound for the allowable solution variables. (DP(N)) If the algorithm chooses X(I)

	- C               : Vector that controls the step length adjustment.  The suggested value for all elements is 2.0. (DP(N))

	- IPRINT          : controls printing inside SA. (INT), 0 - 3, greater levels of printing.

	- ISEED1          : The first seed for the random number generator RANMAR. 0 <= ISEED1 <= 31328. (INT)

	- ISEED2          : The second seed for the random number generator RANMAR. 0 <= ISEED2 <= 30081.

	- T               : On input, the initial temperature. See Goffe et al. for advice. On output, the final temperature. (DP)

	- VM              : The step length vector. On input it should encompass the region of interest given the starting value X.

	- T      : On input, the initial temperature. See Goffe et al. for advice. On output, the final temperature. (DP)

	- VM     : The step length vector. On input it should encompass the region of interest given the starting value X.

	- dat1   : Data object to store simulations.

	- res1   : Dynprob object to store solutions.

	- mn     : Data moments vector.

	- W      : Variance of data moments.

	- momv   : Simulated moment vector.

	- randus : Matrix of random numbers for simulations.

	- des    : String name for saving files.

	- object : Objective function to be minimized.

+ diffevol : Differential evolution code. Finds minimum of objective function, objec, of the form in objec.cpp.

	Arguments:

	- minb        : minimum bounds of problem

	- maxb        : maximum bounds of problem

	- cr          : crossover probability (differential evolution parameter)

	- f           : differential weight (differential evolution parameter)

	- npp         : number of population members per processor

	- maxgen      : maximum number of generations

	- m           : data moment matrix

	- W           : data variance-covariance matrix

	- seedn       : random seed for reproducability of evolution

	- comm_id     : MPI communicator group

	- vali        : dynamic problem operator to keep track of solutions

	- randus      : random vector for simulations

	- nob         : number of

	- des         : filename string to save file to

	- objec       : objective function

	- dati        : data object to store simulations

	- init_points : initial points matrix of size n x p where p is the parameter vector size and 0 <= n <= npp * ncores

	- init_obj    : matrix size n, 0 <= n <= npp * ncores for objective function values for initial points. If you want the program to recalculate objective functions can be left blank.




+ neldmead : Nelder-Mead code.




## OBJECTIVES

The general format to go with optimization, these take inputs:

- par    : parameter vector.

- m      : data moments/EPF.

- W      : data variance covariance matrix.

- res1   : dynprob object to store reslts of VFI.

- dat1   : data object for simulation.

- randus : random uniform matrix for simulations.

- nrun   : number of simulations.

- mn     : vector to store moments from model.


and output a double to be minimied. I'm working on getting the slope EPFs back, currently working are:

+ objec_epf_bspline : Basis spline EPF

+ objec_epf_quad    :	Polynomial EPF

+ objec_mom         : Moments objective
