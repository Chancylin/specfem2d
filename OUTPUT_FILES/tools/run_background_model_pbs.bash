#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher_solver_ak135f_global
#PBS -j oe
#PBS -o job.o
##PBS -q sandy

###########################################################
# USER PARAMETERS

## 256 CPUs ( 16*16 ), walltime 4 hour
#PBS -l nodes=6:ppn=8,walltime=4:00:00

###########################################################
module purge
module load intel/12.1.3
module load openmpi/1.4.4-intel-v12.1

nproc=48
#warning, if you have processors more than 99, you may 
#want to change the 'cat' processing to get the correct
#format of the file name so it could combine all the files

current_dir=$PWD
echo $current_dir
SEM_dir='/scratch/l/liuqy/lcx/specfem2d_self_couple/specfem2d'
cd $SEM_dir
mesh_dir='./DATA/earth_model_ak135'

#obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid
###########################################################
echo `date`
#echo 'runing ak135f background model'
##do the copying
#cp $mesh_dir/global_homo/Par_file ./DATA/
#cp $mesh_dir/global_homo/SOURCE ./DATA/
#cp $mesh_dir/global_homo/STATIONS ./DATA/

echo 'runing ak135f model with LLSVP'
#do the copying
cp $mesh_dir/global_hete/Par_file ./DATA/
cp $mesh_dir/global_hete/SOURCE ./DATA/
cp $mesh_dir/global_hete/STATIONS ./DATA/

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '6,6c\F' ./DATA/switch_solver
sed -i '10,10c\F' ./DATA/switch_solver
sed -i '14,14c\F' ./DATA/switch_solver
sed -i '18,18c\F' ./DATA/switch_solver

#run mesh
echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_ak135f_hete


#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_ak135f_hete

echo `data`
echo 'finish ak135f background running'



