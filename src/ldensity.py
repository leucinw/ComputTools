
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import subprocess
import numpy as np

def main():
  usage = ''' Usage: python ldensity.py xxx.arc '''
  arcfile = sys.argv[1]
  if len(sys.argv) != 2:
    sys.exit(usage)
  if not os.path.isfile(arcfile):
    sys.exit(f"{arcfile} does not exist!")
  
  massdict = {"H": 1.007941, "C":12.01074, "N":14.006703, "O":15.999405, "S":32.0648, "P":30.973761998, \
              "F": 18.99840, "CL":35.4529, "BR":79.9035, "NA":22.989769, "K":39.0983, "MG":24.3051, "ZN":65.38}
  
  natom = subprocess.check_output(["head", "-n1", arcfile]).split()[0]
  headstr = "head -n%s %s > tmp.xyz "%(int(natom)+2, arcfile)
  os.system(headstr)
  grepstr = "grep ' 90.000000' %s > lattice.dat"%(arcfile)
  os.system(grepstr)
  atoms = np.loadtxt("tmp.xyz", usecols=(1,), unpack=True, dtype='str', skiprows=2)
  a,b,c = np.loadtxt("lattice.dat", usecols=(0,1,2), unpack=True, dtype='float')
  os.system("rm -rf tmp.xyz lattice.dat")
  mass = 0.0
  navogadro = 6.02214086 #x10**23
  for atom in atoms:
    mass+=massdict[atom[:1].upper()]
  sumdens = 0.0
  volume = (a*b*c).mean()
  dens = 1000*10.0*(mass/navogadro)/volume
  print("%20s%10.1f kg/m^3" %("Average Density:", dens))
  return

if __name__ == "__main__":
  main()
