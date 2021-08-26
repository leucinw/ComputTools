
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# find three eigenvalues from polarizability tensor 
# calculated from Gaussian output "Exact polarizability"
# or Psi4 output "Dipole Polarizability" tensor 

import argparse
import numpy as np

def diag(output,filetype):
  bohr2angstrom = 0.529177249
  psi4 = False
  x1 = []; x2 = []; x3 = []
  lines=open(output).readlines()
  for line in lines:
    if filetype == 0: #Gaussian
      if "Exact polarizability" in line:
        data = line.split()
        x1 = [float(data[-6]), float(data[-5]), float(data[-3])]
        x2 = [float(data[-5]), float(data[-4]), float(data[-2])]
        x3 = [float(data[-3]), float(data[-2]), float(data[-1])]
    elif filetype == 1: #Psi4
      if "Dipole Polarizability" in line:
        idx = lines.index(line) + 7 
        psi4 = True
  if psi4:
    d = lines[idx+0].split()
    x1 = [float(d[1]), float(d[2]), float(d[3])]
    d = lines[idx+1].split()
    x2 = [float(d[1]), float(d[2]), float(d[3])]
    d = lines[idx+2].split()
    x3 = [float(d[1]), float(d[2]), float(d[3])]

  RED = '\033[91m'
  ENDC = '\033[0m'
  if (x1 != []) and (x2 != []) and (x3 != []):
    a = np.array([x1,x2,x3])
    a = a * bohr2angstrom**3.0
    eigvals = np.linalg.eigvals(a)
    eigvals.sort()
    print("%30s%10.4f%10.4f%10.4f"%(output, eigvals[0], eigvals[1], eigvals[2]))
  else:
    print(RED + "%30s   Polarizability tensor not found"%output + ENDC)
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'input',  nargs='+', help = "Input files: gaussian output files of polarizability calculation", required=True)  
  parser.add_argument('-f', dest = 'filetype', help = "0: Gaussian; 1:Psi4", default=0, type=int)  
  args = vars(parser.parse_args())
  inps = args["input"]
  ftype = args["filetype"]
  print("%30s%10s%10s%10s" %("Filename", "a1(A.^3)", "a2(A.^3)", "a3(A.^3)"))
  [diag(inp,ftype) for inp in inps]
  return

if __name__ == "__main__":
  main()
