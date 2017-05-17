#!/bin/bash
#SBATCH --nodes=50
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=14
#SBATCH --time=4:00:00

#echo commands to stdout
set -x


#move to working directory
cd /pylon1/se4s82p/rkahn/EPF/C

#copy input files to LOCAL file storage
#srun -N $SLURM_NNODES --ntasks-per-node=1 \
#  sh -c 'cp /pylon2/se4s82p/rkahn/input.${SLURM_PROCID} $LOCAL'

cp /home/rkahn/EPF/C/monte_epfq2 /pylon1/se4s82p/rkahn/EPF/C

#run MPI program
mpirun -print-rank-map -n $SLURM_NTASKS -genv \
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -genv I_MPI_PIN_DOMAIN=omp ./monte_epfq2

#copy output files to persistent space
