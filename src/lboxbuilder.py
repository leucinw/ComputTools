
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import argparse
import os
import numpy as np

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz',dest = 'xyz', required=True)  
  parser.add_argument('-key',dest = 'key', required=True)  
  parser.add_argument('-size',dest = 'size', required=True)  
  parser.add_argument('-density',dest = 'density', required=True)  
  parser.add_argument('-tinkerpath',dest = 'path', required=True)  
  args = vars(parser.parse_args())

  xyz = args["xyz"]
  key = args["key"]
  size = float(args["size"])
  path = args["path"]
  density = float(args["density"])

  xyzedit = os.path.join(path, "xyzedit.x")
  if not os.path.isfile(xyzedit):
    xyzedit = os.path.join(path, "xyzedit")
   
  massdict = {"H": 1.007941, "C":12.01074, "N":14.006703, "O":15.999405, "S":32.0648, "P":30.973761998, \
              "F":18.998403163, "CL":35.4529, "BR":79.9035,"NA":22.98976928, "K":39.0983, "MG":24.3051, "ZN":65.38 }
  
  navogadro = 6.02214086 #x10**23
  
  mass = 0.0
  atoms = np.loadtxt(xyz, usecols=(1,), skiprows=1, dtype="str", unpack=True)
  for atom in atoms:
    mass += massdict[atom.upper()]
  volume = size**3.0
  nmolecule = (int(navogadro*volume/(mass/density)*0.1))
  with open("tmp.in", "w") as f:
    f.write("19\n%s\n%s %s %s\nY\n"%(nmolecule, size,size,size))
  cmdstr = "%s %s -k %s < tmp.in "%(xyzedit, xyz, key)
  os.system(cmdstr)
  if not os.path.isfile("%s_2"%xyz):
    print("Box builder failed for %s !"%xyz)
  else:
    os.system("mv %s_2 %s-%sA.xyz"%(xyz, xyz[:-4], int(size)))
    
if __name__ == "__main__":
  main()
