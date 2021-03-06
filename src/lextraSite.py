import numpy as np
import sys,os

'''Usage: generate lone pair "ATOMs" for anisotropic van der Waals interaction
Author: Chengwen Liu
Date: 09/30/2017 (first version)
      05/07/2020 (second version)'''

'''Designed for the following atoms:
   1. sp3 Oxygen, e.g. H2O, CH3OH, CH3-O-CH3
   2. sp2 Oxygen, e.g. (CH3)2C=O
   3. sp3 Nitrogen, e.g. NH3
   4. sp2 Nitrogen, e.g. N=C
   5. sp3 Sulphur, e.g. R-SH, R-S-R, DMSO 
''' 

LPAtomList = ["O", "N", "S"]
typeLP =   {"O":998, "N":999, "S":997}
symbolLP = {"O":"ho", "N":'hn', "S":'hs'}
angleList = {"O":135.0,"N":130.0, "S":130.0} #Only for sp3 O,N,S with 2 connections
lengthList = {"O":1.00, "N":1.00, "S":1.00}  #distance ONS_LP

def main(txyz, atomidx, hybrad):
  connections = []
  atoms = np.loadtxt(txyz, usecols=(1,), skiprows=1, dtype='str', unpack=True) 
  x,y,z = np.loadtxt(txyz, usecols=(2,3,4), skiprows=1, dtype='float', unpack=True)
  txyzlines = open(txyz).readlines()
  for line in txyzlines[1:]:
    connect = []
    dd = line.split()
    for d in dd[6:]:
      connect.append(int(d))
    connections.append(connect) 

  atomidx = int(atomidx)-1
  coordONS = np.array([x[atomidx], y[atomidx], z[atomidx]])
  '''sp3 hybridization, connected 3 atoms''' #NH3,MeNH2,AcNH2
  if (int(hybrad) == 3) and (len(connections[atomidx])==3):
    atomidx1 = connections[atomidx][0]-1
    atomidx2 = connections[atomidx][1]-1
    atomidx3 = connections[atomidx][2]-1
    coordConnect = np.array([[x[atomidx1], y[atomidx1], z[atomidx1]], \
                             [x[atomidx2], y[atomidx2], z[atomidx2]], \
                             [x[atomidx3], y[atomidx3], z[atomidx3]],])
    H_ONS_1 = coordONS - coordConnect[0]
    H_ONS_2 = coordONS - coordConnect[1]
    H_ONS_3 = coordONS - coordConnect[2]

    #Unit vector of trisector
    trisect = H_ONS_1/np.linalg.norm(H_ONS_1) + H_ONS_2/np.linalg.norm(H_ONS_2) + H_ONS_3/np.linalg.norm(H_ONS_3)
    trisect = trisect/np.linalg.norm(trisect)

    #Lone pair position vectors
    length = lengthList[atoms[atomidx]]
    point1 =coordONS + length*trisect 
    iatom = len(atoms)+1
    with open(txyz+"_2", "w") as f:
      f.write("%3s Generated by lextraSite.py\n"%(len(atoms)+1))
      for n in range(1,len(txyzlines)):
        if (n-1) == atomidx:
          f.write(txyzlines[n][:-1] + "  %4s\n"%(len(atoms)+1))
        else:
          f.write(txyzlines[n])
      f.write("%6s%4s %12.6f%12.6f%12.6f%6s  %4s\n"%(iatom,symbolLP[atoms[atomidx]],point1[0],point1[1],point1[2],typeLP[atoms[atomidx]],atomidx+1))
  return 

main(sys.argv[1], sys.argv[2], sys.argv[3])

