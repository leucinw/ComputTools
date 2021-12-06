
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys

def pdb2xyz():
  lines = open(pdb).readlines() 
  allresiduelib = bioresiduelib + solvresiduelib
  for line in lines: 
    if ("TER" not in line) and ("CRYST" not in line) and ("END" not in line):
      resname = line[17:20]
      if resname not in allresiduelib:
        supported = False
        sys.exit(f"{resname} in {pdb} is not supported")
  
  # biomolecules
  chain = 'A B C D E F G'.split() 
  pointer = 0
  with open("bio.pdb", 'w') as f:
    for line in lines: 
      if ("TER" in line):
        pointer += 1
      if ("TER" not in line) and ("CRYST" not in line) and ("END" not in line):
        resname = line[17:20]
        if resname in bioresiduelib: 
          line = line[:21] + chain[pointer] + line[22:]
          f.write(line)
  chains = []
  biolines = open("bio.pdb").readlines()
  for line in biolines:
    if line[21] not in chains:
      chains.append(line[21])
  nchains = len(chains)
  if nchains > 1:
    cmd = "rm bio.xyz* bio.seq*; $TINKER8C8/pdbxyz.x bio.pdb ALL amoebabio18.prm"
  else:
    cmd = "rm bio.xyz* bio.seq*; $TINKER8C8/pdbxyz.x bio.pdb amoebabio18.prm"
  os.system(cmd)
  # ions and water
  if os.path.isfile("bio.xyz"):
    bioxyzlines = open("bio.xyz").readlines()
    index = int(bioxyzlines[0].split()[0])
    with open("sol.xyz", 'w') as f:
      for line in lines:
        if ("TER" not in line) and ("CRYST" not in line) and ("END" not in line):
          resname = line[17:20]
          atomname = line[13:16]
          x = float(line[30:38])
          y = float(line[38:46])
          z = float(line[46:54])
          if resname in solvresiduelib:
            if resname == "Na+":
              atomtype = "352"
              index += 1
              connections = []
              atom = "Na"
            if resname in ["K+ ", "POT"]:
              atomtype = "353"
              index += 1
              connections = []
              atom = "K"
            if resname in ["Cl-", "CLA"]:
              atomtype = "361"
              index += 1
              connections = []
              atom = "Cl"
            if (resname in ["WAT", "TIP"]) and (atomname in ["O ", "OH2"]):
              atomtype = "349"
              index += 1
              connections = [str(index + 1), str(index + 2)]
              atom = "O"
            if (resname in ["WAT", "TIP"]) and (atomname == "H1 "):
              atomtype = "350"
              index += 1
              connections = [str(index - 1)]
              atom = "H"
            if (resname in ["WAT", "TIP"]) and (atomname == "H2 "):
              atomtype = "350"
              index += 1
              connections = [str(index - 2)]
              atom = "H"
            formatted = f"{index:>6d} {atom:2s}{x:12.6f}{y:12.6f}{z:12.6f} {atomtype:5s}" + '   '.join(connections) + "\n"
            f.write(formatted)
    if "CRYST" in lines[0]:
      with open("bio.xyz", 'w') as f:
        d = lines[0].split()
        f.write(f"{index:>5d} Generated by pdb2txyz.py\n")
        f.write(f"{float(d[1]):12.6f}{float(d[2]):12.6f}{float(d[3]):12.6f}{float(d[4]):12.6f}{float(d[5]):12.6f}{float(d[6]):12.6f}\n")
        for line in bioxyzlines[1:]:
          f.write(line)
  lines = open('bio.xyz').readlines()[1:]
  with open("bio.xyz_1", 'w') as f:
    f.write(f"{index:>5d} Generated by pdb2txyz.py\n")
    for line in lines:
      f.write(line)
  catcmd = f"cat bio.xyz_1 sol.xyz > {pdb.replace('pdb', 'xyz')}; rm bio.xyz_1"
  os.system(catcmd)
  return

if __name__ == "__main__":
  global pdb
  pdb = sys.argv[1]
  global bioresiduelib, solvresiduelib
  amberbiolib = ["DA ", "DC ", "DT ", "DG ", "DA5", "DC5", "DT5", "DG5", "DA3", "DC3", "DT3", "DG3"]
  charmmbiolib = ["GUA", "CYT", "THY", "ADE"]
  solvresiduelib = ["WAT", "Na+", "Cl-", " K+", "POT", "CLA", "TIP"]
  bioresiduelib = amberbiolib + charmmbiolib
  pdb2xyz()