
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import numpy as np

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\033[93m'

def integerizer():
  atomtypes, elements = np.loadtxt(xyz, usecols=(5,1), unpack=True, skiprows=1, dtype='str')
  type_element_dict = {}
  for t, e in zip(atomtypes,elements):
    if t not in type_element_dict.keys():
      type_element_dict[t] = e
  type_charge_dict = {}
  for line in open(prm).readlines():
    d = line.split()
    if ("multipole " in line) and (d[1] in atomtypes) and ("#" not in line):
      if d[1] not in type_charge_dict:
        type_charge_dict[d[1]] = (float(d[-1]))
      else:
        if type_charge_dict[d[1]] != (float(d[-1])):
          print(YELLOW + f"Warning: two sets of multipoles differ for {d[1]}" + ENDC)
  chgsum = 0 
  for a in atomtypes:
    chgsum += type_charge_dict[a]
 
  heavy_atom_degeneracy = {}
  for k in type_charge_dict.keys():
    if type_element_dict[k].upper()[0] != 'H':
      heavy_atom_degeneracy[k] = list(atomtypes).count(k)
  
  reschg = chgsum - chg 
  
  probchg = 0.00001
  if reschg > 0:
    minus = True
  else:
    minus = False

  print(GREEN + f" Orginal Charge: {chgsum:10.5f}" + ENDC)
  print(GREEN + f"  Target Charge: {chg:10.5f}" + ENDC)
  
  if minus:
    while (reschg > 0.00001):
      for a in heavy_atom_degeneracy.keys():
        if (reschg < 0.00001):
          break
        deg = heavy_atom_degeneracy[a]
        trychg = deg*probchg
        reschg = reschg - trychg

        if reschg < 0:
          reschg += trychg
          pass
        else:
          type_charge_dict[a] -= probchg
  else:
    while (reschg < -0.00001):
      for a in heavy_atom_degeneracy.keys():
        if (reschg > -0.00001):
          break
        deg = heavy_atom_degeneracy[a]
        trychg = deg*probchg
        reschg = reschg + trychg

        if reschg > 0:
          reschg -= trychg
          pass
        else:
          type_charge_dict[a] += probchg
    
  reschg = float("%8.5f"%reschg)
  for a in heavy_atom_degeneracy.keys():
    deg = heavy_atom_degeneracy[a]
    if deg == 1:
      type_charge_dict[a] -= reschg
      break

  chgsum = 0.0 
  for a in atomtypes:
    chgsum += type_charge_dict[a]
  print(GREEN + "Total charge integerized to: %8.5f electrons"%chgsum + ENDC)
  return type_charge_dict

def integerizerH():
  atomtypes, elements = np.loadtxt(xyz, usecols=(5,1), unpack=True, skiprows=1, dtype='str')
  type_element_dict = {}
  for t, e in zip(atomtypes,elements):
    if t not in type_element_dict.keys():
      type_element_dict[t] = e
  type_charge_dict = {}
  for line in open(prm).readlines():
    d = line.split()
    if ("multipole " in line) and (d[1] in atomtypes) and ("#" not in line):
      if d[1] not in type_charge_dict:
        type_charge_dict[d[1]] = (float(d[-1]))
      else:
        if type_charge_dict[d[1]] != (float(d[-1])):
          print(YELLOW + f"Warning: two sets of multipoles differ for {d[1]}" + ENDC)
  chgsum = 0 
  for a in atomtypes:
    chgsum += type_charge_dict[a]
  hydrogen_atom_degeneracy = {}
  for k in type_charge_dict.keys():
    if type_element_dict[k].upper()[0] == 'H':
      hydrogen_atom_degeneracy[k] = list(atomtypes).count(k)
  
  reschg = chgsum - chg 
  
  probchg = 0.00001
  if reschg > 0:
    minus = True
  else:
    minus = False

  print(GREEN + f" Orginal Charge: {chgsum:10.5f}" + ENDC)
  print(GREEN + f"  Target Charge: {chg:10.5f}" + ENDC)
  
  if minus:
    while (reschg > 0.00001):
      for a in hydrogen_atom_degeneracy.keys():
        if (reschg < 0.00001):
          break
        deg = hydrogen_atom_degeneracy[a]
        trychg = deg*probchg
        reschg = reschg - trychg

        if reschg < 0:
          reschg += trychg
          pass
        else:
          type_charge_dict[a] -= probchg
  else:
    while (reschg < -0.00001):
      for a in hydrogen_atom_degeneracy.keys():
        if (reschg > -0.00001):
          break
        deg = hydrogen_atom_degeneracy[a]
        trychg = deg*probchg
        reschg = reschg + trychg

        if reschg > 0:
          reschg -= trychg
          pass
        else:
          type_charge_dict[a] += probchg
    
  reschg = float("%8.5f"%reschg)
  for a in hydrogen_atom_degeneracy.keys():
    deg = hydrogen_atom_degeneracy[a]
    if deg == 1:
      type_charge_dict[a] -= reschg
      break

  chgsum = 0 
  for a in atomtypes:
    chgsum += type_charge_dict[a]
  print(GREEN + "Total charge integerized to: %8.5f electrons"%chgsum + ENDC)
  return type_charge_dict
def writePRM(type_charge_dict):
  lines = open(prm).readlines()
  with open(prm + "_1", 'w') as f:
    for line in lines:
      d = line.split()
      if ("multipole " in line) and (d[1] in type_charge_dict) and ("#" not in line):
        if "%10.5f"%float(d[-1]) != "%10.5f"%type_charge_dict[d[1]]:
          line = '     '.join(d[:-1] + ["%10.5f"%type_charge_dict[d[1]]]) + "\n"
          f.write(line)
        else:
          f.write(line)
      else:
        f.write(line)
  return 
  
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz', dest = 'xyz', required=True, help="tinker xyz")
  parser.add_argument('-prm', dest = 'prm', required=True, help="tinker key")
  parser.add_argument('-chg', dest = 'chg', required=True, help="total charge", type=float)
  parser.add_argument('-hyd', dest = 'hyd', help="Hydrogen only")
  args = vars(parser.parse_args())
  global xyz,prm,chg
  xyz = args["xyz"]
  prm = args["prm"]
  chg = args["chg"]
  hyd = args["hyd"]
  charges = integerizer()
  if hyd.upper() == "Y":
    charges = integerizer()
  else:
    charges = integerizerH()
  writePRM(charges)
  return

if __name__ == "__main__":
  main()