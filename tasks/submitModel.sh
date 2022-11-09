#!/bin/bash
#BSUB -o log.out
#BSUB -e log.err
#BSUB -q $2
#BSUB -R "span[ptile=36]"
#BSUB -n $1
#BSUB -K

export I_MPI_PLATFORM="skx"
export I_MPI_EXTRA_FILE_SYSTEM=1
export I_MPI_HYDRA_BOOTSTRAP="lsf"
export I_MPI_HYDRA_BRANCH_COUNT=$(( $( echo "${LSB_MCPU_HOSTS}" | wc -w ) / 2 ))
export I_MPI_HYDRA_COLLECTIVE_LAUNCH=1

mpiexec.hydra $3 
touch $4/shell_complete

