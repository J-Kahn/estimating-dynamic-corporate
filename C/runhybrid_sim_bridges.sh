#!/bin/bash
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 14

#echo commands to stdout
set -x


#move to working directory
cd /pylon5/se4s82p/rkahn/estimating-dynamic-corporate/C

export OMP_NUM_THREADS=14

cp /home/rkahn/estimating-dynamic-corporate/C/monte_sim /pylon5/se4s82p/rkahn/estimating-dynamic-corporate/C

./monte_sim
