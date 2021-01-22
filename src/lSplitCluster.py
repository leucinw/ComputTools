

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

"""split a tinker xyz dimer cluster file to its monomer xyzs"""

import os
import sys
import argparse

def split(pool):
  mono1 = [pool[0]]
  pool.remove(pool[0])
  for m in mono1:
    dd = lines[m].split()[6:]
    for d in dd:
      if (int(d) not in mono1):
        mono1.append(int(d))
        if int(d) in pool:
          pool.remove(int(d))
  return sorted(mono1),pool

def writeMono(splitted):
  for i in range(len(splitted)):
    fname = prefix + "_m%02d.xyz"%(i+1)
    with open(fname, "w") as f:
      f.write("%3s\n"%len(splitted[i]))
      for j in range(len(splitted[i])):
        dd = lines[splitted[i][j]].split()
        delta = max(splitted[i]) - len(splitted[i])
        newline = "%3s%3s%12.6f%12.6f%12.6f%5s %s \n"%(j+1, dd[1], float(dd[2]), float(dd[3]), float(dd[4]), dd[5], "  ".join([str(int(d)-delta) for d in dd[6:]]))  
        f.write(newline)
  return

def main():
  global txyz
  global lines
  global prefix 
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'inputs',  nargs='+', help = "Input files",    required=True)  
  args = vars(parser.parse_args())
  inps = args["inputs"]
  for txyz in inps:
    prefix, _ = os.path.splitext(txyz) 
    lines = open(txyz).readlines()
    natom = len(lines)
    pool = []
    for i in range(1, natom, 1):
      pool.append(i)
    splitted = []
    while pool != []:
      r = split(pool)
      mono1, pool = r[0], r[1]
      splitted.append(mono1) 
    writeMono(splitted)
  return

if len(sys.argv) == 1:
  print('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
else:
  main()
