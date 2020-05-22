
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import numpy as np
import subprocess
import sys

''' Usage: python ldensity arcfile'''

massdict = {"H": 1.007941, "C":12.01074, "N":14.006703, "O":15.999405, "S":32.0648, "P":30.973761998, \
            "F":18.998403163, "CL":35.4529, "BR":79.9035,"NA":22.98976928, "K":39.0983, "MG":24.3051, "ZN":65.38 }
arcfile = sys.argv[1]
natom = subprocess.check_output(["head", "-n1", arcfile]).split()[0]
headstr = "head -n%s %s > tmp.xyz "%(int(natom)+2, arcfile)
subprocess.run(headstr,shell=True)
grepstr = "grep ' 90.000000' %s > lattice.dat"%(arcfile)
subprocess.run(grepstr,shell=True)
atoms = np.loadtxt("tmp.xyz", usecols=(1,), unpack=True, dtype='str', skiprows=2)
lattices = np.loadtxt("lattice.dat", usecols=(0,), unpack=True, dtype='float')
mass = 0.0
navogadro = 6.02214086 #x10**23
for atom in atoms:
  mass+=massdict[atom.upper()]
for i in range(len(lattices)):
  volume = lattices[i]**3
  dens = 10.0*(mass/navogadro)/volume
  print("%20s%10.4f"%("ARC_Frame_%s"%i, dens))
