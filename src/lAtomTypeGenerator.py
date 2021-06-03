
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import argparse
import numpy as np
import os
import sys

# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

def genAtomType(txyz, mode):
  fname, _ = os.path.splitext(txyz)
  lines = open(txyz).readlines()
  atomnums, elements, atomtyps = np.loadtxt(txyz, usecols=(0,1,5,), unpack=True, dtype="str", skiprows=1)
  num_ele_dict = dict(zip(atomnums, elements))
  typ_ele_dict = {}
  for typ, ele in zip(atomtyps, elements):
    if typ not in typ_ele_dict:
      typ_ele_dict[typ] = ele
  lines = open(txyz).readlines()[1:]
  connect = []
  for line in lines:
    dd = line.split()
    connect.append([d for d in dd[6:]])
  with open(f"{fname}.type", 'w') as f:
    for line in lines:
      dd = line.split()
      typ = dd[5]
      if mode == "0":
        typ_str = typ_ele_dict[typ]+"*"
      elif mode == "1": 
        strs = []
        for s in dd[6:]:
          strs.append(num_ele_dict[s])

        if (typ_ele_dict[typ].upper() == 'H'):
          nb = len(connect[int(dd[6])-1])
          typ_str = typ_ele_dict[typ] + ''.join(strs) + str(nb)
        else:
          typ_str = typ_ele_dict[typ] + ''.join(sorted(strs))
      else:
        sys.exit(RED + "run with -h to see the usage" + ENDC) 
      f.write(f"{dd[0]:3s} {dd[5]:5s} {typ_str:10s}\n")
      #print(GREEN + f"{dd[0]:3s} {dd[5]:5s} {typ_str:10s}" + ENDC)
  print(GREEN + f"{txyz}" + ENDC)
  return

 
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'txyz', help = "tinker xyz file", required=True)  
  parser.add_argument('-m', dest = 'mode', help = "0: generic type; 1: type based on the first neighboring atoms;", default = '1')  
  args = vars(parser.parse_args())
  txyz = args["txyz"]
  mode = args["mode"]
  genAtomType(txyz, mode)
  return

if __name__ == "__main__":
  main()
