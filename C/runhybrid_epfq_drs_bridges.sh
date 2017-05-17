#!/bin/bash
#SBATCH --nodes=25
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=14
#SBATCH --time=12:00:00

#echo commands to stdout
set -x


#move to working directory

export OMP_NUM_THREADS=14

#copy input files to LOCAL file storage
#srun -N $SLURM_NNODES --ntasks-per-node=1 \


#run MPI program
mpirun -print-rank-map -n $SLURM_NTASKS -genv \
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -genv I_MPI_PIN_DOMAIN=omp ./monte_epfq_drs

#copy output files to persistent space
