
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

from readChemFile import *
from geomFunction import *

import sys

usage = "python bondangle.py filename atom1 atom2 [atom3]"

if (len(sys.argv) !=4) and (len(sys.argv) !=5):
  print(usage)
else:
  if ".txyz" in sys.argv[1]:
    data = readTXYZ(sys.argv[1])
  if ".xyz" in sys.argv[1]:
    data = readXYZ(sys.argv[1])
  # bond
  if len(sys.argv) == 4:
    idx_1 = int(sys.argv[2])
    idx_2 = int(sys.argv[3])
    coord_1 = data[1][idx_1-1]
    coord_2 = data[1][idx_2-1]
    print(distance(coord_1, coord_2))

  # angle
  else:  
    idx_1 = int(sys.argv[2])
    idx_2 = int(sys.argv[3])
    idx_3 = int(sys.argv[4])
    coord_1 = data[1][idx_1-1]
    coord_2 = data[1][idx_2-1]
    coord_3 = data[1][idx_3-1]
    print(angle(coord_1, coord_2, coord_3))
