#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher_solver_ak135f_step4
#PBS -j oe
#PBS -o job.o
#PBS -q sandy

###########################################################
# USER PARAMETERS

## 256 CPUs ( 16*16 ), walltime 4 hour
#PBS -l nodes=3:ppn=16,walltime=8:00:00

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

#Step 4
echo 'runing step 4'
#do the copying
cp $mesh_dir/global_homo/Par_file ./DATA/
cp $mesh_dir/global_homo/SOURCE ./DATA/
cp $mesh_dir/global_homo/STATIONS ./DATA/

#change the nproc in Par_file. Make sure line 19 is for the nproc
awk -F'=' '{if (NR==19) {OFS="="; $2=nproc}; print}' nproc=" $nproc"  ./DATA/Par_file > ./DATA/temp_par
mv ./DATA/temp_par ./DATA/Par_file

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '4,4c\F' ./DATA/switch_solver
sed -i '8,8c\F' ./DATA/switch_solver
sed -i '12,12c\F' ./DATA/switch_solver
sed -i '16,16c\T' ./DATA/switch_solver

echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step4

#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step4

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./global_hete_reconst
