
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import argparse,sys
import numpy as np

# color
RED = '\33[91m'
GREEN = '\33[92m'
ENDC = '\033[0m'

def main():
  if len(sys.argv) == 1:
    sys.exit(RED + "run with -h option to see the usage"+ ENDC)

  parser = argparse.ArgumentParser()
  parser.add_argument('probe', default="charge", type=str.lower)
  parser.add_argument('-i', dest = 'txyz', required=True, help="tinker xyz file not containing the probe ion")  
  parser.add_argument('-p', dest = "prob_atom", nargs='+', required=True, type=int, help="probe atom followed by directional atoms, e.g. 1 2 3: use bisector of 1-2 and 1-3") 
  parser.add_argument('-d', dest = 'distance', nargs='+', required=True, type=float, help="distance from the probing atom, e.g. 3 4 5: three distances at 3, 4 and 5 A")
  args = vars(parser.parse_args())
  txyz = args["txyz"]
  patom = args["prob_atom"]
  dists  = args["distance"] 
  probe  = args["probe"] 

  # read atoms and coordinate from txyz file
  def readTXYZ(txyz):
    atoms  = np.loadtxt(txyz, usecols=(1,), dtype='str', unpack=True, skiprows=1)
    coords = np.loadtxt(txyz, usecols=(2,3,4), dtype='float', unpack=False, skiprows=1)
    return atoms,coords

  # calculate distance between two atoms 
  def distance(coord1, coord2):
    return np.sqrt(np.square(np.array(coord1)-np.array(coord2)).sum()) 

  
  def writeCOM(atoms, coords, fname):
    with open(fname,"w") as f:
      nproc = "%nproc=8\n" 
      memory = "%mem=20GB\n"
      chk = "%chk=" + "%s\n"%fname.replace(".com", ".chk")
      keywords = "#p MP2/Aug-cc-pvtz SP Density=MP2 SCF=Save Charge NoSymm MaxDisk=100GB \n"
      if probe == "ion":
        keywords = "#p MP2/Aug-cc-pvtz SP Density=MP2 SCF=Save NoSymm MaxDisk=100GB \n"
      comment = "SP job with a positive charge\n"
      chgspin = " 0 1  \n"
      if probe == "ion":
        chgspin = " 1 1  \n"
      f.write(chk + memory + nproc + keywords + "\n" + comment + "\n" + chgspin)
      for atom, coord in zip(atoms, coords):
        f.write("%3s %10.6f%10.6f%10.6f\n"%(atom, coord[0], coord[1], coord[2]))
    return
 
  # write COM and XYZ files 
  atoms, coords = readTXYZ(txyz)
  vector = np.zeros(3)
  lines = open(txyz).readlines()
  for j in range(1,len(patom)):
    vector += (coords[patom[0]-1] - coords[patom[j]-1])
  vector = vector/np.linalg.norm(vector)
  for dist in dists:
    v = coords[patom[0]-1] + dist*vector
    fname = txyz.split(".")[0] + "_%02dA.com"%dist
    writeCOM(atoms, coords, fname)
    with open(fname, "a") as f:
      if probe == "charge":
        f.write("\n")
        f.write(" %10.6f%10.6f%10.6f%4.1f\n\n"%(v[0], v[1], v[2], 1.0))
      elif probe == "ion":
        f.write("%3s %10.6f%10.6f%10.6f\n"%("Na", v[0], v[1], v[2]))
        f.write("\n")
      else:
        sys.exit(RED + "supported probes: charge or ion" + ENDC)
      print(GREEN + "Generated %s"%fname + ENDC)

    with open(fname.replace(".com", ".xyz"), "w") as f:
      f.write("  %4s\n"%(len(atoms)+1))
      for line in lines[1:]:
        f.write(line)
      if probe == "charge":
        f.write("  %4s%3s  %12.6f%12.6f%12.6f      999\n"%(len(atoms)+1, "X", v[0], v[1], v[2]))
      elif probe == "ion":
        f.write("  %4s%3s  %12.6f%12.6f%12.6f        7\n"%(len(atoms)+1, "Na", v[0], v[1], v[2]))
      else:
        sys.exit(RED + "supported probes: charge or ion" + ENDC)
      print(GREEN + "Generated %s"%fname.replace(".com", ".xyz") + ENDC)
  return 

if __name__ == "__main__":
  main()
