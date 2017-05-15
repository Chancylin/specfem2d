#!/bin/bash

nproc=2

current_dir=$PWD
echo $current_dir

SEM_dir='/home/lcx/specfem2d_devep/specfem2d'
cd $SEM_dir

mesh_dir='./DATA/quick_test'

#Step 4
echo 'runing step 4'
#do the copying
cp $mesh_dir/global/Par_file ./DATA/
cp $mesh_dir/global/SOURCE ./DATA/
cp $mesh_dir/global/STATIONS ./DATA/

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '4,4c\F' ./DATA/switch_solver
sed -i '8,8c\F' ./DATA/switch_solver
sed -i '12,12c\F' ./DATA/switch_solver
sed -i '16,16c\T' ./DATA/switch_solver

echo 'do the meshing'
mpirun.mpich -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step4

#run solver
echo 'do the solver'
mpirun.mpich -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step4

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./global_hete_reconst
