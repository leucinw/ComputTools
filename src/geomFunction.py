
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

def dihedral(coord1, coord2, coord3, coord4):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array(coord1) 
    p1 = np.array(coord2) 
    p2 = np.array(coord3) 
    p3 = np.array(coord4) 

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x)) 
