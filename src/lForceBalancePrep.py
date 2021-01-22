
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

'''Prepare an interaction.txt that can be used as a ForceBalance interaction energy target'''

import os
import sys
import string
import numpy as np

inputfile = sys.argv[1]
dimerxyzs = np.loadtxt(inputfile, usecols=(0), dtype="str", unpack=True, skiprows=1) 
energies,weights = np.loadtxt(inputfile, usecols=(1,2), dtype="float", unpack=True, skiprows=1) 

alphabet = string.ascii_lowercase
ndimer = len(dimerxyzs)
if (3*ndimer <= 26):
  shortNameList = list(alphabet)
elif (26 < 3*ndimer <= 26*26):
  for i in list(alphabet):
    for j in list(alphabet):
      shortNameList.append(i+j)
elif (26*26 < 3*ndimer <= 26*26*26):
  for i in list(alphabet):
    for j in list(alphabet):
      for k in list(alphabet):
        shortNameList.append(i+j+k)
else:
  print("Too many xyz file provided?")

allxyzs = []
for xyz in dimerxyzs:
  if "_1.00.xyz" in xyz:
    allxyzs.append(xyz.replace(".xyz", "_m01.xyz"))
    allxyzs.append(xyz.replace(".xyz", "_m02.xyz"))

BEnames = ["BE_" + x for x in shortNameList[len(allxyzs): len(allxyzs) + ndimer]]

allxyzs = allxyzs + list(dimerxyzs)
ntotal = len(allxyzs)
systems = shortNameList[0:ntotal] 
xyz_system = dict(zip(allxyzs, systems))

with open("interactions.txt",'w') as f:
  f.write('$global\nkeyfile interactions.key\nenergy_unit kilocalories_per_mole\n$end\n\n')
  for syst, clst in zip(systems,allxyzs): 
    f.write('$system\nname %s\ngeometry %s\n$end\n\n'%(syst, clst))
  for i in range(len(dimerxyzs)):
    k1 = dimerxyzs[i][:-8] + "1.00_m01.xyz"
    k2 = dimerxyzs[i][:-8] + "1.00_m02.xyz"
    equ = xyz_system[dimerxyzs[i]] + " - "  + xyz_system[k1] + " - " + xyz_system[k2]
    f.write('$interaction\nname %s\nequation %s\nenergy  %s\nweight %s\n$end\n\n'%(BEnames[i], equ, energies[i], weights[i]))
