#!/usr/bin/env python

# translate the second molecule so that we have a certain distance 
# defined by two lists of atoms
# refAtomsFirstMol:
# refAtomsSecondMol

from geomFunction import *

def f_translate(firstMolCoord, secondMolCoord, refAtomsFirstMol, refAtomsSecondMol, distance_frac):
  coord_1 = firstMolCoord 
  coord_2 = secondMolCoord
  refCoord_1 = []
  refCoord_2 = []
  nfirst = len(firstMolCoord)
  for i in refAtomsFirstMol:
    refCoord_1.append(firstMolCoord[i-1])
  for i in refAtomsSecondMol:
    refCoord_2.append(secondMolCoord[i-1-nfirst])

  geocent_1 = geom_center(refCoord_1)
  geocent_2 = geom_center(refCoord_2)

  geocent_dist = distance(geocent_1, geocent_2)
  geocent_vector = [geocent_2[0]- geocent_1[0], 
                    geocent_2[1]- geocent_1[1], 
                    geocent_2[2]- geocent_1[2],]
  trans_vector = [0, 0, 0]
  # this is the speciall case where two molecules are the same
  if geocent_dist == 0: 
    for i in range(len(coord_2)):
      coord_2[i][0] += 1.0
    geocent_2 = geom_center(coord_2)
    geocent_dist = distance(geocent_1, geocent_2)
    geocent_vector = [geocent_2[0]- geocent_1[0], 
                      geocent_2[1]- geocent_1[1], 
                      geocent_2[2]- geocent_1[2],]

  #dist_ratio = target_dist/geocent_dist
  dist_ratio = distance_frac 
  trans_vector = [geocent_vector[0]*(dist_ratio-1.0),
                  geocent_vector[1]*(dist_ratio-1.0),
                  geocent_vector[2]*(dist_ratio-1.0), ]
  # translate second molecule
  for i in range(len(coord_2)):
    for j in range(3):
      coord_2[i][j] += trans_vector[j] 

  totCoord = coord_1 + coord_2

  # return the translated coordinates
  return totCoord
 
