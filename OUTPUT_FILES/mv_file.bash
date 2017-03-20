#!/bin/bash
mkdir $1
echo 'move simulation output to ' $1
cp ../DATA/Par_file ../DATA/SOURCE ../DATA/STATIONS $1 
mv image00* test2.rec* wavefield* Ux_file_single.bin Uz_file_single.bin Database* source.txt $1
echo 'finish'
