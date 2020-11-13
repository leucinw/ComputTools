
'''Prepare an interaction.txt that can be used as a ForceBalance interaction energy target'''

import os,sys

clusters = []
monomers = []
energies = []
systems  = []
BEnames  = []
weights  = []

lines = open(sys.argv[1]).readlines()
for line in lines:
  dd = line.split()
  clusters.append(dd[0])
  monomers.append(dd[2:-2])
  energies.append(dd[-2])
  weights.append(dd[-1])
  systems.append(dd[1])
  BEnames.append("BE_"+dd[0].split(".xyz")[0])

mono_sys_dict = { "Methane.xyz"   : "Meth",
                  "Ethane.xyz"    : "Eth",
                  "Propane.xyz"   : "Prop",
                  "Methanol.xyz"  : "MeOH",
                  "Ethanol.xyz"   : "EtOH",
                  "1-Propanol.xyz": "nPropOH",
                  "2-Propanol.xyz": "isoPropOH",
                  "dimer01_mono01.xyz" : "dimer01_mono01",
                  "dimer02_mono01.xyz" : "dimer02_mono01",
                  "dimer03_mono01.xyz" : "dimer03_mono01",
                  "dimer04_mono01.xyz" : "dimer04_mono01",
                  "dimer05_mono01.xyz" : "dimer05_mono01",
                  "dimer06_mono01.xyz" : "dimer06_mono01",
                  "dimer07_mono01.xyz" : "dimer07_mono01",
                  "dimer08_mono01.xyz" : "dimer08_mono01",
                  "dimer09_mono01.xyz" : "dimer09_mono01",
                  "dimer10_mono01.xyz" : "dimer10_mono01",
                  "dimer11_mono01.xyz" : "dimer11_mono01",
                  "dimer12_mono01.xyz" : "dimer12_mono01",
                  "dimer13_mono01.xyz" : "dimer13_mono01",
                  "dimer14_mono01.xyz" : "dimer14_mono01",
                  "dimer15_mono01.xyz" : "dimer15_mono01",
                  "dimer16_mono01.xyz" : "dimer16_mono01",
                  "dimer17_mono01.xyz" : "dimer17_mono01",
                  "dimer18_mono01.xyz" : "dimer18_mono01",
                  "dimer19_mono01.xyz" : "dimer19_mono01",
                  "dimer20_mono01.xyz" : "dimer20_mono01",
                  "dimer01_mono02.xyz" : "dimer01_mono02",
                  "dimer02_mono02.xyz" : "dimer02_mono02",
                  "dimer03_mono02.xyz" : "dimer03_mono02",
                  "dimer04_mono02.xyz" : "dimer04_mono02",
                  "dimer05_mono02.xyz" : "dimer05_mono02",
                  "dimer06_mono02.xyz" : "dimer06_mono02",
                  "dimer07_mono02.xyz" : "dimer07_mono02",
                  "dimer08_mono02.xyz" : "dimer08_mono02",
                  "dimer09_mono02.xyz" : "dimer09_mono02",
                  "dimer10_mono02.xyz" : "dimer10_mono02",
                  "dimer11_mono02.xyz" : "dimer11_mono02",
                  "dimer12_mono02.xyz" : "dimer12_mono02",
                  "dimer13_mono02.xyz" : "dimer13_mono02",
                  "dimer14_mono02.xyz" : "dimer14_mono02",
                  "dimer15_mono02.xyz" : "dimer15_mono02",
                  "dimer16_mono02.xyz" : "dimer16_mono02",
                  "dimer17_mono02.xyz" : "dimer17_mono02",
                  "dimer18_mono02.xyz" : "dimer18_mono02",
                  "dimer19_mono02.xyz" : "dimer19_mono02",
                  "dimer20_mono02.xyz" : "dimer20_mono02",
                }

with open("interactions.txt",'w') as f:
  f.write('$global\nkeyfile interactions.key\nenergy_unit kilocalories_per_mole\n$end\n\n')
  for syst, clst in zip(systems,clusters): 
    f.write('$system\nname %s\ngeometry %s\n$end\n\n'%(syst, clst))
  for mono in mono_sys_dict:
    if os.path.isfile(mono): 
      f.write('$system\nname %s\ngeometry %s\n$end\n\n'%(mono_sys_dict[mono], mono))
  for i in range(len(systems)):
    equ = systems[i]
    for j in range(len(monomers[i])):
      equ += " - " + mono_sys_dict[monomers[i][j]] + " "
    f.write('$interaction\nname %s\nequation %s\nenergy  %s\nweight %s\n$end\n\n'%(BEnames[i], equ, energies[i], weights[i]))

