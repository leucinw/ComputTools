
#===================================
# 	     Chengwen Liu              #
# 	   liuchw2010@gmail.com        #
# 	University of Texas at Austin  #
#===================================

# find three eigenvalues from polarizability tensor 
# calculated from Gaussian output "Exact polarizability"

# usage: polarSolver.py Gau.out

import numpy
from numpy import *
from sympy import *
import sys

bohr2angstrom = 0.529177249
x1=[];x2=[];x3=[]

def eigenvalue(filename):
  lines=open(filename).readlines()
  for line in lines:
    if "Exact polarizability" in line:
      data=line.split()
      x1=[float(data[-6]), float(data[-5]), float(data[-3])]
      x2=[float(data[-5]), float(data[-4]), float(data[-2])]
      x3=[float(data[-3]), float(data[-2]), float(data[-1])]
  
  a = Matrix([x1,x2,x3])
  a = a * bohr2angstrom**3.0
  b = a.diagonalize()
  print("%30s%10.4f%10.4f%10.4f"%(filename,float(str(b[-1][0]).split()[0]),float(str(b[-1][4]).split()[0]),float(str(b[-1][8]).split()[0])))


eigenvalue(sys.argv[1])
