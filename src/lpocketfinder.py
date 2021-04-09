
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import sys
import os
import numpy as np

## gold.conf: regular gold.conf file 
## sub.sh: should be /path/to/your/gold/gold_auto gold.conf

pro_pdb = sys.argv[1] 
lig_mol = sys.argv[2] 

usage = "python lpocketfinder.py protein.pdb ligand.mol2"

def grid():
  x = []
  y = []
  z = []
  lines = open(pro_pdb).readlines()
  for line in lines:
    dd = line.split()
    ## coordinate line
    if (len(dd) == 11) and (dd[0] == "ATOM"):
      x.append(float(dd[5]))
      y.append(float(dd[6]))
      z.append(float(dd[7]))
  x = np.array(x)
  y = np.array(y)
  z = np.array(z)

  x_min, x_max = int(10*round(x.min()/10.0)), int(10*round(x.max()/10.0)) 
  y_min, y_max = int(10*round(y.min()/10.0)), int(10*round(y.max()/10.0)) 
  z_min, z_max = int(10*round(z.min()/10.0)), int(10*round(z.max()/10.0))
  for i in range(x_min, x_max+1, 20):
    for j in range(y_min, y_max+1, 20):
      for k in range(z_min, z_max+1, 20):
        submit(i,j,k)
  return

def submit(x,y,z):
  if not os.path.isfile("gold.conf"):
    sys.exit("please provide a gold.conf file")
  with open("temp.conf","w") as f:
    for line in open("gold.conf").readlines():
      if "origin" in line:
        line = "origin = %5.1f %5.1f %5.1f \n"%(float(x), float(y), float(z))
      if "protein_datafile" in line:
        line = f"protein_datafile = ../{pro_pdb}\n"
      if "ligand_data_file" in line:
        line = f"ligand_data_file ../{lig_mol} 10 \n"
      f.write(line)
  dirname = "grid_%03dx%03dx%03d"%(x,y,z)
  os.system(f"mkdir {dirname}")
  os.system(f"cp temp.conf {dirname}/gold.conf")
  os.system(f"cp sub.sh ./{dirname}")
  workdir = os.path.join(currdir, dirname)
  os.system(f"cd {workdir} && ./sub.sh")
  return

global currdir
currdir = os.getcwd()
grid()
