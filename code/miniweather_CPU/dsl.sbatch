#!/bin/bash
#SBATCH -J DSL
#SBATCH -o DSL.out
#SBATCH -e DSL.err
#SBATCH -p cnCPU
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --ntasks-per-socket=8
#SBATCH --exclusive

module purge;
module load compilers/intel/psxe_2020_u4.all basekits/zlib/1.2.11-psxe_2020_u4 basekits/szip/2.1.1-psxe_2020_u4 basekits/curl/7.79.1-psxe_2020_u4 basekits/hdf5/1.12.1-psxe_2020_u4_mpi basekits/netcdf_c/4.8.1-psxe_2020_u4_mpi basekits/netcdf_f/4.5.3-psxe_2020_u4_mpi compilers/cmake/3.21.4

export FI_PROVIDER="^ofi_rxm"
ulimit -s unlimited
ulimit -c unlimited

#srun hostname | sort -u > hosts.%j

mpirun ./app
