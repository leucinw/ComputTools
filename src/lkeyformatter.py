
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import sys
import os
import argparse

'''write a nice key/prm file from several input key/prm files; 
   formated key/prm reduces the possibilities of making mistakes '''

# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', dest = 'forcefield', help = "Force Field header", choices=["AMOEBA", "AMOEBA+", "NONE"], default="AMOEBA+")  
  parser.add_argument('-p', dest = 'potentials', nargs='+', help = "Potential term(s) to write", default = [])  
  parser.add_argument('-k', dest = 'tinkerkeys', nargs='+', help = "tinker .key or .prm file(s)",default = [])  
  args = vars(parser.parse_args())
  ff   = args["forcefield"]
  pots = args["potentials"]
  keys = args["tinkerkeys"]
  keywords = ['atom','bond','angle','strbnd','opbend','torsion','vdw','polarize','chgpen','chgtrn','bndcflux','angcflux','multipole']
  if pots != []:
    newkeywords = []
    for pot in pots:
      if not pot.lower() in keywords:
        sys.exit(RED + "Potential term is not recognized!" + ENDC)
      else:
        newkeywords.append(pot.lower())
    keywords = newkeywords
  if keys == []:
    sys.exit(RED + "Please provide input key(s)! Or run with -h to see the usage"+ ENDC)
  else:
    for key in keys:
      if not os.path.isfile(key):
        sys.exit(RED + f"{key} you provided does not exist!" + ENDC)

  if ff == "AMOEBA":
    header = "/home/liuchw/poltype-latest/ParameterFiles/amoebabio18_header.prm"
    os.system(f"cat {header}")
  elif ff == "AMOEBA+":
    header = "/home/liuchw/amoebaplus.data/converged.prm/header.prm"
    os.system(f"cat {header}")
  else:
    pass

  for keyword in keywords:
    if keyword != 'multipole':
      grepstr = "grep '%s ' %s --no-filename"%(keyword, '  '.join([key for key in keys]))
      os.system(grepstr)
  
  if 'multipole' in keywords:
    for key in keys:
      lines = open(key).readlines()
      for i in range(len(lines)):
        d = lines[i].split()
        if (len(d) > 1) and (d[0] == 'multipole'):
          print(lines[i+0][:-1])
          print(lines[i+1][:-1])
          print(lines[i+2][:-1])
          print(lines[i+3][:-1])
          print(lines[i+4][:-1])

if __name__ == "__main__":
  main()
