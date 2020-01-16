#!/bin/bash

# To run on other cards, GTX1070 etc. if you decide to use this, comment out the above blocks 
source /home/liuchw/.openmm.amoeba.cuda8.0 
export filename=k_water
export CUDA_VISIBLE_DEVICES=0 # device number
export OPENMMHOME=/home/liuchw/Softwares/tinkers/Tinker-OpenMM/AMOEBA/Tinker8.5/source
T=298.15

# NVT simulations, 2fs time step, save every 2 ps, 0.1ns in total, 298 K
nohup $OPENMMHOME/dynamic_omm.x $filename.xyz -k $filename.key 50000 2.0 2.0 2 $T N > $filename.log 2>err.log & 
