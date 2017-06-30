#!/bin/bash

nproc=2

current_dir=$PWD
echo $current_dir

SEM_dir='/home/lcx/specfem2d_devep/specfem2d'
cd $SEM_dir

mesh_dir='./DATA/conceptual_model'

#Step 2
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
sed -i '6,6c\T' ./DATA/switch_solver
sed -i '10,10c\F' ./DATA/switch_solver
sed -i '14,14c\F' ./DATA/switch_solver
sed -i '18,18c\F' ./DATA/switch_solver

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
