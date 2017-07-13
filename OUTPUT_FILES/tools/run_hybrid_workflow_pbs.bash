#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher_solver_earth_like_hybrid_psv_offset_test
#PBS -j oe
#PBS -o job.o
##PBS -q sandy

###########################################################
# USER PARAMETERS

## 256 CPUs ( 16*16 ), walltime 4 hour
#PBS -l nodes=6:ppn=8,walltime=8:00:00

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
#Step 1
echo `date`
echo 'runing step 1'
#do the copying
cp $mesh_dir/local_hete/Par_file ./DATA/
cp $mesh_dir/local_hete/SOURCE ./DATA/
cp $mesh_dir/local_hete/STATIONS ./DATA/

#change the nproc in Par_file. Make sure line 19 is for the nproc
awk -F'=' '{if (NR==19) {OFS="="; $2=nproc}; print}' nproc=" $nproc"  ./DATA/Par_file > ./DATA/temp_par 
mv ./DATA/temp_par ./DATA/Par_file

rm ./DATA/boundary_points*

#set the switch
sed -i '2,2c\T' ./DATA/switch_solver
sed -i '6,6c\F' ./DATA/switch_solver
sed -i '10,10c\F' ./DATA/switch_solver
sed -i '14,14c\F' ./DATA/switch_solver
sed -i '18,18c\F' ./DATA/switch_solver

#run mesh
echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step1

#run solver to get the boundary_points

#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step1

#move the result into the directory


#combine all the boundary_points prifiles into one
num=`expr $nproc - 1`
file_total='./DATA/boundary_points_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    if [ $i -lt 10 ];then
       file_proc="./DATA/boundary_points0000"$i
       cat $file_proc >> $file_total
    elif [ $i -lt 100 ];then
       file_proc="./DATA/boundary_points000"$i
       cat $file_proc >> $file_total
    fi
done

#######################################################################################

#Step 2
echo `date`
echo 'runing step 2'
#do the copying
cp $mesh_dir/global_homo/Par_file ./DATA/
cp $mesh_dir/global_homo/SOURCE ./DATA/
cp $mesh_dir/global_homo/STATIONS ./DATA/

#make directory to store the bd_info
mkdir -p ./OUTPUT_FILES/bg_record/
mkdir -p ./OUTPUT_FILES/bg_record/

#change the nproc in Par_file. Make sure line 19 is for the nproc
awk -F'=' '{if (NR==19) {OFS="="; $2=nproc}; print}' nproc=" $nproc"  ./DATA/Par_file > ./DATA/temp_par 
mv ./DATA/temp_par ./DATA/Par_file

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '6,6c\T' ./DATA/switch_solver
sed -i '10,10c\F' ./DATA/switch_solver
sed -i '14,14c\F' ./DATA/switch_solver
sed -i '18,18c\F' ./DATA/switch_solver

echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step2

#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step2

##################################################################
#combine the profiles from the global model into one
num=`expr $nproc - 1`
file_total='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    if [ $i -lt 10 ];then
       file_proc="./OUTPUT_FILES/bg_record/elastic_pnts_profile0000"$i
       cat $file_proc >> $file_total
    elif [ $i -lt 100 ];then
       file_proc="./OUTPUT_FILES/bg_record/elastic_pnts_profile000"$i
       cat $file_proc >> $file_total
    fi
done

file_total='./OUTPUT_FILES/bg_record/acoustic_pnts_profile_total'

rm $file_total
for i in $(seq 0 $num)
do
    echo $i
    if [ $i -lt 10 ];then
       file_proc="./OUTPUT_FILES/bg_record/acoustic_pnts_profile0000"$i
       cat $file_proc >> $file_total
    elif [ $i -lt 100 ];then
       file_proc="./OUTPUT_FILES/bg_record/acoustic_pnts_profile000"$i
       cat $file_proc >> $file_total
    fi
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

#change the nproc in Par_file. Make sure line 19 is for the nproc
awk -F'=' '{if (NR==19) {OFS="="; $2=nproc}; print}' nproc=" $nproc"  ./DATA/Par_file > ./DATA/temp_par 
mv ./DATA/temp_par ./DATA/Par_file

#make directory to store the bd_info
mkdir -p ./OUTPUT_FILES/reconst_record/
mkdir -p ./OUTPUT_FILES/reconst_record/

#set the switch
sed -i '2,2c\F' ./DATA/switch_solver
sed -i '6,6c\F' ./DATA/switch_solver
sed -i '10,10c\T' ./DATA/switch_solver
sed -i '14,14c\T' ./DATA/switch_solver
sed -i '18,18c\F' ./DATA/switch_solver

echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step3

#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step3
##################################################################
#combine the profiles from the global model into one
num=`expr $nproc - 1`
file_total='./OUTPUT_FILES/reconst_record/elastic_pnts_profile_total'
rm $file_total

for i in $(seq 0 $num)
do
    echo $i
    if [ $i -lt 10 ];then
       file_proc="./OUTPUT_FILES/reconst_record/elastic_pnts_profile0000"$i
       cat $file_proc >> $file_total
    elif [ $i -lt 100 ];then
       file_proc="./OUTPUT_FILES/reconst_record/elastic_pnts_profile000"$i
       cat $file_proc >> $file_total
    fi
done

file_total='./OUTPUT_FILES/reconst_record/acoustic_pnts_profile_total'

rm $file_total
for i in $(seq 0 $num)
do
    echo $i
    if [ $i -lt 10 ];then
       file_proc="./OUTPUT_FILES/reconst_record/acoustic_pnts_profile0000"$i
       cat $file_proc >> $file_total
    elif [ $i -lt 100 ];then
       file_proc="./OUTPUT_FILES/reconst_record/acoustic_pnts_profile000"$i
       cat $file_proc >> $file_total
    fi
done
##################################################################

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./local_hete

#######################################################################################

#Step 4

cd $SEM_dir

echo `date`
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
sed -i '6,6c\F' ./DATA/switch_solver
sed -i '10,10c\F' ./DATA/switch_solver
sed -i '14,14c\F' ./DATA/switch_solver
sed -i '18,18c\T' ./DATA/switch_solver

echo 'do the meshing'
mpirun -np $nproc ./bin/xmeshfem2D > ./OUTPUT_FILES/output_mesher_step4

#run solver
echo 'do the solver'
mpirun -np $nproc ./bin/xspecfem2D > ./OUTPUT_FILES/output_solver_step4

#move the result and the configuration files
cd ./OUTPUT_FILES
./mv_file.bash ./global_hete_reconst


echo `date`
echo "all is done, cheers"
