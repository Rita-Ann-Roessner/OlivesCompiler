#!/bin/bash

#SBATCH --time=20:00:00 
#SBATCH --nodes=2
#SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4
#SBATCH --account=ehm@v100
#SBATCH --output=out.%j 
#SBATCH --error=err.%j 
#SBATCH --job-name=test 
#SBATCH --exclusive
#SBATCH -C mps           # turn on CUDA MPS to enable multiple processes per GPU
#SBATCH --hint=nomultithread    # 1 MPI process per physical core (no hyperthreading)
#SBATCH --qos=qos_gpu-t3

module purge
module load gcc/12.2.0
module load cuda/12.4.1
module load openmpi/4.1.5
#module load gromacs/2023.4-cuda
module load gromacs/2024.2-cuda-plumed

source ~/.bashrc
conda activate cgcompiler-env

srun --cpu-bind=verbose -n 20 ./gpu_bind_multi4.sh  python -u run.py -n ${SLURM_CPUS_PER_TASK}
