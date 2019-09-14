#!/usr/bin/env python

# distance between two atoms
# input the two coordinates[x,y,z]
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
