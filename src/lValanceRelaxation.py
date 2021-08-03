
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import numpy as np

def prepare(xyz):
  cmdstr = f"python {converterdir}/lToTXYZ.py -i {xyz} -t xyz"
  os.system(cmdstr)
  with open("tinker.key", 'w') as f:
    f.write("parameters ttt.prm\n")
    f.write("bondterm only\n")
    f.write("openmp-threads 8\n")
  return

def genprm():
  elements, types = np.loadtxt(txyz, usecols=(1,5), dtype='str', unpack=True, skiprows=1) 
  headers = open(f"{amoebaheaddir}/amoebabio18_header.prm").readlines()
  tmp = []
  with open("ttt.prm",'w') as f:
    for line in headers:
      f.write(line)
    for e,t in zip(elements, types):
      if t not in tmp:
        f.write(f'atom    {t:>6s} {t:>6s} {e:>6s} "{e}" {ptable[e][0]:>5d} {ptable[e][1]:10.3f} {ptable[e][2]:10d}\n')
        tmp.append(t)
  cmdstr = f"{tinkerdir}/analyze.x {txyz} -k tinker.key > ana_err.dat"
  os.system(cmdstr)
  lines = open("ana_err.dat").readlines()
  bondlist = [] 
  for line in lines:
    if "Bond   " in line:
      d = line.split()
      a1,a2 = int(d[-2]), int(d[-1])
      if a1 > a2:
        a1, a2 = a2, a1
      if [a1,a2] not in bondlist:
        bondlist.append([a1,a2])
  with open("ttt.prm", 'a') as f:
    f.write("# ASSIGN\n")
    for b in bondlist:
      f.write(f"bond {b[0]:>6d} {b[1]:>6d} 0.0 0.0 \n")
  cmdstr = f"python {assignprmdir}/lAssignAMOEBAplusPRM.py -xyz {txyz} -key ttt.prm -potent bonded -konly no"
  os.system(cmdstr)
  newkeyname = f"{txyz.replace('.txyz', '_new.key')}"
  bonds = []
  lines = open("ttt.prm").readlines()
  for l in lines:
    if "bond " not in l:
      bonds.append(l)
  lines = open(newkeyname)
  for l in lines:
    if "bond " in l:
      bonds.append(l)
  with open("ttt.prm", 'w') as f:
    for b in bonds:
      f.write(b)
  return

def optimize():
  cmdstr = f"{tinkerdir}/minimize.x {txyz} -k tinker.key 0.5"
  os.system(cmdstr)
  return

def main():
  global txyz
  global ptable
  global tinkerdir
  global converterdir
  global amoebaheaddir
  global assignprmdir

  tinkerdir="/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8"
  converterdir="/home/liuchw/bin"
  amoebaheaddir="/home/liuchw/shared/forcefield"
  assignprmdir="/home/liuchw/bin/AssignPRM"

  ptable = {'H':[1,   1.008, 1], 'C':[6,  12.011, 4], 'N':[7, 14.007, 3],  'O':[8,  15.999, 2], \
            'S':[16, 32.066, 2], 'P':[15, 30.974, 4], 'F':[9, 18.998, 1], 'Cl':[17, 35.453, 1]}
  xyz = sys.argv[1]
  txyz = xyz.replace('xyz', 'txyz')
  if os.path.isfile(txyz):
    print(f"Warning: using the existing txyz file {txyz} to assign parameters")
  else:
    prepare(xyz)
  genprm()
  optimize()

if __name__ == "__main__":
  main()