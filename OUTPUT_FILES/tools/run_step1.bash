#!/bin/bash

nproc=2

current_dir=$PWD
echo $current_dir
SEM_dir='/home/lcx/specfem2d_devep/specfem2d'
cd $SEM_dir
mesh_dir='./DATA/quick_test'

#Step 1
echo 'runing step 1'
#do the copying
cp $mesh_dir/local_hete/Par_file ./DATA/
cp $mesh_dir/local_hete/SOURCE ./DATA/
cp $mesh_dir/local_hete/STATIONS ./DATA/

rm ./DATA/boundary_points*

#set the switch
sed -i '2,2c\T' ./DATA/switch_solver
sed -i '4,4c\F' ./DATA/switch_solver
sed -i '8,8c\F' ./DATA/switch_solver
sed -i '12,12c\F' ./DATA/switch_solver
sed -i '16,16c\F' ./DATA/switch_solver

#run mesh
echo 'do the meshing'
mpirun.mpich -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step1

#run solver to get the boundary_points

#run solver
echo 'do the solver'
mpirun.mpich -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step1

#move the result into the directory


#combine all the boundary_points prifiles into one
num=`expr $nproc - 1`
file_total='./DATA/boundary_points_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    file_proc="./DATA/boundary_points0000"$i
    cat $file_proc >> $file_total
done

