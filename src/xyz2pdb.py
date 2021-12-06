
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import argparse
import subprocess
import numpy as np
import concurrent.futures 
from datetime import datetime

def readxyz():
  lines = open(xyz).readlines()
  natoms = int(lines[0].split()[0])
  nlines = len(lines)
  delta = nlines - natoms
  if delta == 1:
    ftype = "TXYZ"
  elif delta == 2:
    ncols = len(lines[2].split())
    if ncols == 4:
      ftype = "XYZ"
    else:
      ftype = "TXYZ"
  else:
    sys.exit("Could not guess the file type")
 
  pbc = []
  if ftype == "TXYZ":
    atoms = np.loadtxt(xyz, usecols=(1,), dtype='str', unpack=True, skiprows=delta)  
    xs,ys,zs = np.loadtxt(xyz, usecols=(2,3,4), dtype='float', unpack=True, skiprows=delta)
    if delta == 2:
      pbc = [float(x) for x in lines[1].split()]
  else:
    atoms = np.loadtxt(xyz, usecols=(0,), dtype='str', unpack=True, skiprows=delta)  
    xs,ys,zs = np.loadtxt(xyz, usecols=(1,2,3), dtype='float', unpack=True, skiprows=delta)
  return atoms, xs, ys, zs, pbc

def handlepdb(pdb):
  atoms, xs, ys, zs, pbc = readxyz()
  natoms = len(atoms)
  lines = open(pdb).readlines()
  nlines = 0
  for line in lines:
    if ("TER" not in line) and ("END" not in line) and ("CRYST1" not in line):
      nlines += 1
  if natoms != nlines:
    sys.exit(f"Unequal # of atoms: {pdb}:{nlines} and {xyz}:{natoms}")
  else:
    with open("result.pdb", 'w') as f:
      idx = 0
      if pbc != []:
        f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (pbc[0], pbc[1], pbc[2], pbc[3], pbc[4], pbc[5]))
      for i in range(len(lines)):
        if ("TER" not in lines[i]) and ("END" not in lines[i]) and ("CRYST1" not in lines[i]):
          d = lines[i]
          atom = d[-3:].strip()
          x = float(d[30:38])
          y = float(d[38:46])
          z = float(d[46:54])
          s = len(atom)
          if atom.lower() != atoms[i-idx].strip()[0:s].lower():
            sys.exit(f"Atom does not match {atom} --- {atoms[i-idx][0:s]}")
          else:
            newline = d[:30] + "%8.3f%8.3f%8.3f"%(xs[i-idx], ys[i-idx], zs[i-idx]) + d[54:] 
            f.write(newline)
        if ("TER" in lines[i]) or ("END" in lines[i]):
          f.write(lines[i])
          idx += 1
  return

if __name__ == "__main__":
  global xyz
  xyz = sys.argv[1]
  pdb = sys.argv[2]
  ## Tinker recognizes certain names
  ions_in_xyzs = {"POT":"K  ", "CLA":"Cl "}
  for ion, ion_ in ions_in_xyzs.items():
    cmd = f"sed 's/{ion}/{ion_}/g' -i {xyz}"
    os.system(cmd)
  ress_in_pdbs = {"CYT":"DC ", "GUA":"DG ", "ADE":"DA ", "THY":"DT "} 
  for res, res_ in ress_in_pdbs.items():
    cmd = f"sed 's/{res}/{res_}/g' -i {pdb}"
    os.system(cmd)
  handlepdb(pdb) 
