#!/bin/bash

# --- SBATCH File Settings ---
#SBATCH --job-name=DiscreteState
#SBATCH --output=output.log
#SBATCH --error=error.log

# --- SBATCH Partition Settings ---
#SBATCH --partition=ffa
#SBATCH --time=03:59:00

# --- SBATCH Computation settings Settings ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=256M


echo "=========================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Working directory: $(pwd)" 
echo "=========================================================="

echo "================ HARDWARE DETALS ========================="
lscpu
echo "=========================================================="

# Modules loading
echo "Loading modules..."
module load oneapi/tbb oneapi/compiler-rt oneapi/umf oneapi/mkl/latest oneapi/compiler-intel-llvm/latest

# OMP Settings
echo "OMP Number of threads: $SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_MAX_ACTIVE_LEVELS=1

# Running computation
echo "Running computation $HOME/DiscreteState/discrete_state"
#$HOME/DiscreteState/debug_discrete_state
$HOME/DiscreteState/discrete_state

echo "=========================================================="
echo "Computation finished."
echo "Job ending."
echo "=========================================================="