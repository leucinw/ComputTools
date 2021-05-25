
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# find three eigenvalues from polarizability tensor 
# calculated from Gaussian output "Exact polarizability"

import sys
import numpy as np

def diag(gauOut):
  bohr2angstrom = 0.529177249
  x1 = []; x2 = []; x3 = []
  lines=open(gauOut).readlines()
  for line in lines:
    if "Exact polarizability" in line:
      data = line.split()
      x1 = [float(data[-6]), float(data[-5]), float(data[-3])]
      x2 = [float(data[-5]), float(data[-4]), float(data[-2])]
      x3 = [float(data[-3]), float(data[-2]), float(data[-1])]
  if (x1 != []) and (x2 != []) and (x3 != []):
    a = np.array([x1,x2,x3])
    a = a * bohr2angstrom**3.0
    eigvals = np.linalg.eigvals(a)
    eigvals.sort()
    print("%30s%10.4f%10.4f%10.4f"%(gauOut, eigvals[0], eigvals[1], eigvals[2]))
  else:
    print("%30s Not found!"%gauOut)
  return

diag(sys.argv[1])
