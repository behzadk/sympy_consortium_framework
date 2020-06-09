#!/bin/bash -l

module unload compilers mpi
module load compilers/gnu/4.9.2
module load python3/3.6

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/ucbtbdk/Scratch/cpp_consortium_sim/model_sel_doe

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/stage/lib

export INCLUDE=$INCLUDE:/lustre/home/ucbtbdk/software/boost_1_65_1/

export C_INCLUDE_PATH=$C_INCLUDE_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/stage/lib
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/libs
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1


export LIBRARY_PATH=$LIBRARY_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/stage/lib
export LIBRARY_PATH=$LIBRARY_PATH:/lustre/home/ucbtbdk/software/boost_1_65_1/libs

# 1. Set num threads
#$ -pe smp 8
export OMP_NUM_THREADS=$(ppn)
echo "$OMP_NUM_THREADS"


# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=48:0:0

# 3. Request 1 gigabyte of RAM (must be an integer)
#$ -l mem=4G
#$ -l tmpfs=10G

# 5. Set the name of the job.
#$ -N BK_three_spec_table

# 5. Set the name of the job.
echo "moving to BK_manu"
cd $HOME/Scratch/cpp_consortium_sim/BK_manu_three_species

./build_myriad.sh


