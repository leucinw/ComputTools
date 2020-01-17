
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


# distance between two atoms
# input the two coordinates[x,y,z]

import numpy as np

def distance(coord1, coord2):
  dist = 0.0
  for i in range(3):
    dist += (coord1[i] - coord2[i])**2.0
  dist = dist**0.5
  return dist

# geometric center of a molecule
# input the coordinates [[x1,y1,z1], [x2,y2,z2], ...]
def geom_center(coord):
  geocent = [0.0, 0.0, 0.0]
  leng = len(coord)
  for i in range(leng):
    for j in range(3):
      geocent[j] += coord[i][j]/leng
  return geocent

# angle of three atoms 
# input the coordinates [x1,y1,z1], [x2,y2,z2], [x3,y3,z3] 
def angle(coord1, coord2, coord3):
  angl = 0.0
  coord1 = np.array(coord1)
  coord2 = np.array(coord2)
  coord3 = np.array(coord3)
  v21 = coord1-coord2
  v23 = coord3-coord2
  dot = np.dot(v21,v23)
  v21norm = sum((coord1-coord2)**2)**0.5
  v23norm = sum((coord3-coord2)**2)**0.5
  angl = np.arccos(dot/(v21norm*v23norm))
  angl = 180.0*angl/np.pi
  return angl  
