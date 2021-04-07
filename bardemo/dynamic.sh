#!/bin/bash

source /home/liuchw/OpenMM-Nov-2020/openmm.bashrc
export filename=liquid
export CUDA_VISIBLE_DEVICES=0 
export OPENMMHOME=/home/liuchw/OpenMM-Nov-2020/tinker/source/
export T=298.15

nohup $OPENMMHOME/dynamic_omm.x $filename.xyz -k $filename.key 50000 2.0 2.0 2 $T N > $filename.log 2>err.log & 
