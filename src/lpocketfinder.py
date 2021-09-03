
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import numpy as np
import concurrent.futures 

## gold.conf: regular gold.conf file 
## sub.sh: should be /path/to/your/gold/gold_auto gold.conf

""" generate necessary files for each grid point """
def generate(grid):
  [x,y,z] = grid
  dirname = "grid_%03dx%03dx%03d"%(x,y,z)
  if not os.path.isfile("gold.conf"):
    sys.exit("please provide a gold.conf file")
  if not os.path.isdir(dirname):
    with open(f"{dirname}.conf","w") as f:
      for line in open("gold.conf").readlines():
        if "origin" in line:
          line = "origin = %5.1f %5.1f %5.1f \n"%(float(x), float(y), float(z))
        if "protein_datafile" in line:
          line = f"protein_datafile = ../{pro_pdb}\n"
        if "ligand_data_file" in line:
          line = f"ligand_data_file ../{lig_mol} 10 \n"
        f.write(line)
    os.system(f"mkdir {dirname}")
    os.system(f"mv {dirname}.conf {dirname}/gold.conf")
    os.system(f"cp sub.sh ./{dirname}")
  return

""" submit the gold job for one grid """
def sub(dirname):
  os.chdir(os.path.join(currdir, dirname))
  os.system("sh sub.sh")
  return

""" generate grids of the protein """
def gengrids():
  x = []; y = []; z = []
  coords = []
  lines = open(pro_pdb).readlines()
  for line in lines:
    dd = line.split()
    if "ATOM " in line: 
      x.append(float(line[30:38]))
      y.append(float(line[38:46]))
      z.append(float(line[46:54]))
      coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
  x = np.array(x)
  y = np.array(y)
  z = np.array(z)
  x_min, x_max = int(10*round(x.min()/10.0)), int(10*round(x.max()/10.0)) 
  y_min, y_max = int(10*round(y.min()/10.0)), int(10*round(y.max()/10.0)) 
  z_min, z_max = int(10*round(z.min()/10.0)), int(10*round(z.max()/10.0))
  grids = []
  for i in range(x_min, x_max+1, 20):
    for j in range(y_min, y_max+1, 20):
      for k in range(z_min, z_max+1, 20):
        grids.append([i,j,k])
  
  jobs = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(generate, g) for g in grids]
    for f in concurrent.futures.as_completed(results):
      jobs.append(f.result())
  
  if grid_file != None:
    if os.path.isfile(grid_file):
      files = []
      for line in open(grid_file).readlines():
        files.append(line.split("\n")[0])
      jobs = []
      with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(sub, f) for f in files]
        for f in concurrent.futures.as_completed(results):
          jobs.append(f.result())
    else:
      sys.exit(f"grid file -- {grid_file} -- does not exist !")
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-pdb', dest = 'protein',  help = "protein PDB file", required=True)  
  parser.add_argument('-mol2', dest = 'ligand',  help = "ligand MOL2 file", required=True)  
  parser.add_argument('-grid', dest = 'grid',  help = "grid.dat file. Optional")  
  args = vars(parser.parse_args())

  global currdir, pro_pdb, lig_mol, grid_file
  pro_pdb = args["protein"]
  lig_mol = args["ligand"]
  grid_file = args["grid"]
  currdir = os.getcwd()
  gengrids()
  return

if __name__ == "__main__":
  main()
