#!/bin/bash

LOCAL_RANK=${MPI_LOCALRANKID} # mpirun Intel MPI
if [ -z "${LOCAL_RANK}" ]; then LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}; fi # mpirun OpenMPI
if [ -z "${LOCAL_RANK}" ]; then LOCAL_RANK=${SLURM_LOCALID}; fi  # srun


if [ $(($LOCAL_RANK%4)) -eq 0 ]
then
  export CUDA_VISIBLE_DEVICES=0
elif [ $(($LOCAL_RANK%4)) -eq 1 ]
then
  export CUDA_VISIBLE_DEVICES=1
elif [ $(($LOCAL_RANK%4)) -eq 2 ]
then
  export CUDA_VISIBLE_DEVICES=2
elif [ $(($LOCAL_RANK%4)) -eq 3 ]
then
  export CUDA_VISIBLE_DEVICES=3
fi

$@