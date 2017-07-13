#!/bin/bash
mkdir $1
echo 'move simulation output to ' $1
cp ../DATA/Par_file ../DATA/SOURCE ../DATA/STATIONS ../DATA/switch_solver $1 
cp ./tools/job.o $1
mv image00* test1.rec* wavefield* Database* source.txt output_solver_step* output_mesher_step* $1
echo 'finish'
