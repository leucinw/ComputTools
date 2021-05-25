#!/bin/bash

source /home/liuchw/OpenMM-Nov-2020/openmm.bashrc
export filename=liquid
export CUDA_VISIBLE_DEVICES=0 
export nextstep
export OPENMMHOME=/home/liuchw/OpenMM-Nov-2020/tinker/source/
T=298.15

nohup $OPENMMHOME/bar_omm.x 1 $filename.arc $T ../$nextstep/$filename.arc $T  >bar1.log 2>err1.log &&
			$OPENMMHOME/bar_omm.x 2 $filename.bar 1 50 1 1 50 1 > freeEnergy.log 2>err2.log &
