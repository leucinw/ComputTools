
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import subprocess
import argparse
import numpy as np

# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

''' Extract a cluster structure from a liquid box; used with tinker.xyz file '''

def distance(p1, p2):
  p1 = np.array(p1)
  p2 = np.array(p2)
  if not hasPBC:
    dist = np.sqrt(np.square(p1-p2).sum()) 
  else:
    dist = np.sqrt(np.square(p1-p2).sum()) 
  return dist

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i',  dest = 'input', help = "liquid box file in tinker xyz format", required=True)  
  parser.add_argument('-r',  dest = 'radii', type=float, help = "radius from the center of solute, in Angstrom", required=True)  
  parser.add_argument('-n1', dest = 'first', type=int, help = "first atom index of the solute molecule")  
  parser.add_argument('-n2', dest = 'final', type=int, help = "final atom index of the solute molecule")  
  args = vars(parser.parse_args())
  inp = args["input"]
  rad = args["radii"]
  at1 = args["first"]
  at2 = args["final"]
  
  global pbc, hasPBC
  lines = open(inp).readlines()
  nline = len(lines)
  natom = int(lines[0].split()[0])
  delta = nline-natom
  if delta == 2:
    hasPBC = True
    d = lines[1].split()
    box = [float(d[0]), float(d[1]), float(d[2])]
  elif delta == 1:
    hasPBC = False
    box = [0., 0., 0.]
  else:
    sys.exit(RED + f"{inp} is not a tinker file" + ENDC)
  
  xs,ys,zs = np.loadtxt(inp, usecols=(2,3,4), skiprows=delta, unpack=True)
  center = [xs[(at1-1):(at2-1)].mean(), ys[(at1-1):(at2-1)].mean(), zs[(at1-1):(at2-1)].mean()]

  for i in range(delta, nline):
    j = i - delta
    if distance(center, [xs[j], ys[j], zs[j]]) <= rad:
      print(lines[i][:-1])

    
  return


if len(sys.argv) == 1:
  print('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
else:
  main()
