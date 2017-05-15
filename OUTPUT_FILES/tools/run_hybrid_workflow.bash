#!/bin/bash

nproc=2

current_dir=$PWD
echo $current_dir
SEM_dir='/home/lcx/specfem2d_devep/specfem2d'
cd $SEM_dir
mesh_dir='./DATA/conceptual_model'

#Step 1
echo `date`
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


#######################################################################################

#Step 2
echo `date`
echo 'runing step 2'
#do the copying
cp $mesh_dir/global/Par_file ./DATA/
cp $mesh_dir/global/SOURCE ./DATA/
cp $mesh_dir/global/STATIONS ./DATA/

#make directory to store the bd_info
mkdir -p ./OUTPUT_FILES/bg_record/
mkdir -p ./OUTPUT_FILES/bg_record/

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '4,4c\T' ./DATA/switch_solver
sed -i '8,8c\F' ./DATA/switch_solver
sed -i '12,12c\F' ./DATA/switch_solver
sed -i '16,16c\F' ./DATA/switch_solver

echo 'do the meshing'
mpirun.mpich -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step2

#run solver
echo 'do the solver'
mpirun.mpich -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step2

##################################################################
#combine the profiles from the global model into one
num=`expr $nproc - 1`
file_total='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    file_proc="./OUTPUT_FILES/bg_record/elastic_pnts_profile0000"$i
    cat $file_proc >> $file_total
done

file_total='./OUTPUT_FILES/bg_record/acoustic_pnts_profile_total'

rm $file_total
for i in $(seq 0 $num)
do
    echo $i
    file_proc="./OUTPUT_FILES/bg_record/acoustic_pnts_profile0000"$i
    cat $file_proc >> $file_total
done
##################################################################

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./global_homo


#######################################################################################

#Step 3

cd $SEM_dir

echo `date`
echo 'runing step 3'
#do the copying
cp $mesh_dir/local_hete/Par_file ./DATA/
cp $mesh_dir/local_hete/SOURCE ./DATA/
cp $mesh_dir/local_hete/STATIONS ./DATA/

#make directory to store the bd_info
mkdir -p ./OUTPUT_FILES/reconst_record/
mkdir -p ./OUTPUT_FILES/reconst_record/

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '4,4c\F' ./DATA/switch_solver
sed -i '8,8c\T' ./DATA/switch_solver
sed -i '12,12c\T' ./DATA/switch_solver
sed -i '16,16c\F' ./DATA/switch_solver

echo 'do the meshing'
mpirun.mpich -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step3

#run solver
echo 'do the solver'
mpirun.mpich -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step3

##################################################################
#combine the profiles from the global model into one
num=`expr $nproc - 1`
file_total='./OUTPUT_FILES/reconst_record/elastic_pnts_profile_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    file_proc="./OUTPUT_FILES/reconst_record/elastic_pnts_profile0000"$i
    cat $file_proc >> $file_total
done

file_total='./OUTPUT_FILES/reconst_record/acoustic_pnts_profile_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    file_proc="./OUTPUT_FILES/reconst_record/acoustic_pnts_profile0000"$i
    cat $file_proc >> $file_total
done
#################################################################

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./local_hete

#######################################################################################

#Step 4

cd $SEM_dir

echo `date`
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

echo `date`
echo "all is done, cheers"
