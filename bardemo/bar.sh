#!/bin/bash

source /home/liuchw/.openmm.amoeba.cuda8.0 
export filename=k_water
export nextstep
export CUDA_VISIBLE_DEVICES=0 # device number
export OPENMMHOME=/home/liuchw/Softwares/tinkers/Tinker-OpenMM/AMOEBA/Tinker8.5/source
T=298.15

nohup $OPENMMHOME/bar_omm.x 1 $filename.arc $T ../$nextstep/$filename.arc $T  >bar1.log 2>err1.log &&
			$OPENMMHOME/bar_omm.x 2 $filename.bar 1 50 1 1 50 1 > freeEnergy.log 2>err2.log &
