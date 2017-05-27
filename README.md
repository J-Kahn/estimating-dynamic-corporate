The makefile in this directory is organized to run one table or set of figures at a time. Before running, if downloaded from github make sure to run **make directories**. The files are meant to run on the bridges supercomputer, though with modifications they should run on other systems. To run on bridges, begin by copying the directory to pylon5 folder and before running python files to invert weight matrices be sure to run **module load python/2.7.11_gcc** so that appropriate python packages will be loaded. The tables can also be generated on bridges, but the version of python differs - tables were originally created using Python 3.5.2, Anaconda 2.4.1 (64-bit). In general, slurm jobs must complete before the next step can be run.

## To replicate Table 3 and Figure 2
1. **make build_large**: builds necessary C++ code and moves it to where it's needed for slurm scripts.
2. **make run_sim_large**: runs the simulations, submitting a single slurm job.
3. **make run_mc_large**: inverts weight matrices for estimation, then submits four slurm jobs for the four columns of Table 3.
4. **make table_large**: makes latex table and figures from the completed monte carlo simulations under the WR folder corresponding to Table 3 and Figure 2.

## To replicate Table 4 and Figure 3
1. **make build_small**: builds necessary C++ code and moves it to where it's needed for slurm scripts.
2. **make run_sim_small**: runs the simulations, submitting a single slurm job.
3. **make run_mc_small**: inverts weight matrices for estimation, then submits four slurm jobs for the four columns of Table 4.
4. **make table_small**: makes latex table and figures from the completed monte carlo simulations under the WR folder corresponding to Table 4 and Figure 3.

## To replicate Table 5
1. **make build_large**: builds necessary C++ code for large sample size.
2. **make build_small**: builds necessary C++ code for small sample size.
3. **make run_sim_large**: runs the simulations for large sample size, submitting a single slurm job.
4. **make run_sim_small**: runs the simulations for small sample size, submitting a single slurm job.
5. **make run_mc_large**: inverts weight matrices for estimation, then submits four slurm jobs for the four columns of Table 5 with a large sample.
6. **make run_mc_small**: inverts weight matrices for estimation, then submits four slurm jobs for the four columns of Table 5 with a small sample.
7. **make table_size**: makes latex table corresponding to Table 5 as well as Figures 2 and 3.

## To replicate Table 6
1. **make build_robust**: builds necessary C++ code and moves it to where it's needed for slurm scripts.
2. **make run_sim_large**: runs the simulations, submitting a single slurm job.
3. **make run_mc_robust**: inverts weight matrices for estimation, then submits three slurm jobs for the three columns of Table 6.
4. **make table_robust**: makes latex table for Table 6.

## To replicate Table 7
1. **make build_fixed**: builds necessary C++ code and moves it to where it's needed for slurm scripts.
2. **make run_sim_fixed**: runs the simulations, submitting a single slurm job.
3. **make run_mc_fixed**: inverts weight matrices for estimation, then submits two slurm jobs for the two columns of Table 7.
4. **make table_fixed**: makes latex table for Table 7.

## To replicate Table 8
1. **make build_drs**: builds necessary C++ code and moves it to where it's needed for slurm scripts.
2. **make run_sim_drs**: runs the simulations, submitting a single slurm job.
3. **make run_mc_drs**: inverts weight matrices for estimation, then submits two slurm jobs for the two columns of Table 8.
4. **make table_drs**: makes latex table for Table 8.

## To replicate Figures 4 and 5
1. **make build_large**: builds necessary C++ code and moves it to where it's needed for slurm scripts for correctly specified model.
2. **make build_mis**: builds necessary C++ code and moves it to where it's needed for slurm scripts for misspecified model.
3. **make run_sim**: runs the simulations, submitting a single slurm job for the correctly specified model.
4. **make run_sim_mis**: runs the simulations, submitting a single slurm job for the misspecified model.
5. **make run_mc_correct**: inverts weight matrices for estimation, then submits two slurm jobs for correctly specified models.
6. **make run_mc_mis**: inverts weight matrices for estimation, then submits two slurm jobs for misspecified models.
7. **make table_mis**: makes Figures 4 and 5 in the WR directory.
