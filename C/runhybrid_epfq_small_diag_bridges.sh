#!/bin/bash
#SBATCH --nodes=25
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=14
#SBATCH --time=12:00:00

#echo commands to stdout
set -x


#move to working directory
cd /pylon5/se4s82p/rkahn/estimating-dynamic-corporate/C

export OMP_NUM_THREADS=14

#copy input files to LOCAL file storage
#srun -N $SLURM_NNODES --ntasks-per-node=1 \
#  sh -c 'cp /pylon2/se4s82p/rkahn/input.${SLURM_PROCID} $LOCAL'

cp /home/rkahn/estimating-dynamic-corporate/C/monte_epfq_small_diag /pylon5/se4s82p/rkahn/estimating-dynamic-corporate/C

#run MPI program
mpirun -print-rank-map -n $SLURM_NTASKS -genv \
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -genv I_MPI_PIN_DOMAIN=omp ./monte_epfq_small_diag

#copy output files to persistent space
cp /pylon5/se4s82p/rkahn/estimating-dynamic-corporate/Results/Trials/* /home/rkahn/estimating-dynamic-corporate/Results/Trials
