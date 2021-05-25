import numpy as np
from convertLib import *
from pybel import *
import sys,os

'''Usage: generate lone pair "ATOMs" for anisotropic exchange repulsion'''
'''Chengwen Liu'''
'''Sep.30th,2017'''

'''Designed for the following atoms:
   1. sp3 Oxygen, e.g. H2O, CH3OH, CH3-O-CH3
   2. sp2 Oxygen, e.g. (CH3)2C=O
   3. sp3 Nitrogen, e.g. NH3
   4. sp2 Nitrogen, e.g. N=C
   5. sp3 Sulphur, e.g. R-SH, R-S-R, DMSO 
''' 

LPAtomList = ["O", "N", "S"]
typeLP = {"O":800, "N":810, "S":820}
symbolLP = {"O":"ho", "N":'hn', "S":'hs'}
angleList = {"O":135.0,"N":130.0, "S":130.0} #Only for sp3 O,N,S
lengthList = {"O":1.00, "N":1.00, "S":1.00} #distance ONS_LP

def main(txyz):
  connections = []
  atoms=[]
  filename = txyz.split(".txyz")[0]
  lines = file(txyz).readlines()
  for line in lines[1:]: 
    data = line.split()
    if data[1] in LPAtomList:  
      connections.append([data[0]] + data[6:])
      atoms.append(data[1]) 
  TXYZ2TXYZ(txyz)
  TXYZ2XYZ(txyz)
  mol = readfile("xyz", filename+".xyz").next()
  iatom = len(lines)
 
  temp = open("tmp","w") 
  for i in xrange(len(connections)):
    coordConnect = []
    hyb = mol.atoms[int(connections[i][0])-1].hyb
    coordONS = mol.atoms[int(connections[i][0])-1].coords
    for j in xrange(1, len(connections[i])):
      coord = mol.atoms[int(connections[i][int(j)])-1].coords
      coordConnect.append(coord)
    coordConnect=np.array(coordConnect)
    ####################################################
    '''sp3 hybridization, connected 2 atoms''' #H2O,HOCH3
    if (hyb == 3) and (len(connections[i])-1==2):
      H_ONS_1 = coordONS - coordConnect[0]
      H_ONS_2 = coordONS - coordConnect[1]

      #Unit vector of bisector
      bisect = H_ONS_1/np.linalg.norm(H_ONS_1) + H_ONS_2/np.linalg.norm(H_ONS_2)
      bisect = bisect/np.linalg.norm(bisect)

      #Unit vector of perpendicular vector
      perpvect = np.cross(H_ONS_1, H_ONS_2)
      perpvect = perpvect/np.linalg.norm(perpvect)

      #Lone pair position vectors
      angle = angleList[atoms[i]]
      length = lengthList[atoms[i]]
      point1 =coordONS + length*np.cos(np.pi*(angle/180.0)/2.0)*bisect + length*np.sin(np.pi*(angle/180.0)/2.0)*perpvect
      point2 =coordONS + length*np.cos(np.pi*(angle/180.0)/2.0)*bisect - length*np.sin(np.pi*(angle/180.0)/2.0)*perpvect

      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point1[0],point1[1],point1[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point2[0],point2[1],point2[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
    ############################################################################################################################################
    '''sp3 hybridization, connected 3 atoms''' #NH3,NH2Ac
    if (hyb == 3) and (len(connections[i])-1==3):
      H_ONS_1 = coordONS - coordConnect[0]
      H_ONS_2 = coordONS - coordConnect[1]
      H_ONS_3 = coordONS - coordConnect[2]

      #Unit vector of trisector
      trisect = H_ONS_1/np.linalg.norm(H_ONS_1) + H_ONS_2/np.linalg.norm(H_ONS_2) + H_ONS_3/np.linalg.norm(H_ONS_3)
      trisect = trisect/np.linalg.norm(trisect)

      #Lone pair position vectors
      length = lengthList[atoms[i]]
      point1 =coordONS + length*trisect 
      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point1[0],point1[1],point1[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
    ############################################################################################################################################
    '''sp2 hybridization, connected 2 atoms'''
    if (hyb == 2) and (len(connections[i])-1==2):
      H_ONS_1 = coordONS - coordConnect[0]
      H_ONS_2 = coordONS - coordConnect[1]

      #Unit vector of bisector
      bisect = H_ONS_1/np.linalg.norm(H_ONS_1) + H_ONS_2/np.linalg.norm(H_ONS_2)
      bisect = bisect/np.linalg.norm(bisect)

      #Lone pair position vectors
      length = lengthList[atoms[i]]
      point1 =coordONS + length*bisect 
      
      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point1[0],point1[1],point1[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
    ############################################################################################################################################
    '''sp2 hybridization, connected 1 atom''' #R1R2C=O
    if (hyb == 2) and (len(connections[i])-1==1):
      betaC = []
      data = lines[ int(connections[i][1])].split()
      for dat in data[6:]:
        if int(dat)!= int(connections[i][0]):
          betaC.append(int(dat))
      alphaC = mol.atoms[int(connections[i][1])-1].coords
      betaC1 = mol.atoms[int(betaC[0])-1].coords
      betaC2 = mol.atoms[int(betaC[1])-1].coords
      
      betaC_alphaC_1 = np.array(alphaC) - np.array(betaC1)
      betaC_alphaC_1 = betaC_alphaC_1/np.linalg.norm(betaC_alphaC_1)

      #Unit vector of beta_alpha vectors
      betaC_alphaC_2 = np.array(alphaC) - np.array(betaC2)
      betaC_alphaC_2 = betaC_alphaC_2/np.linalg.norm(betaC_alphaC_2)
      
      #Lone pair position vectors
      length = lengthList[atoms[i]]
      point1 =coordONS + length*betaC_alphaC_1
      point2 =coordONS + length*betaC_alphaC_2
      
      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point1[0],point1[1],point1[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
      temp.write("%6s%4s%12.6f%12.6f%12.6f%6s%4s\n"%(iatom,symbolLP[atoms[i]],point2[0],point2[1],point2[2],typeLP[atoms[i]],connections[i][0]))
      iatom += 1
    ############################################################################################################################################
    else:
      print(f"No atoms in {txyz} file need to add LP!"
  temp.close()

  of = open("head", "w")
  of.write("%3s Generated by txyzEdit.py\n"%(iatom-1))
  of.close()
  os.system("tail -n%s %s >origtxyz"%((len(lines)-1), txyz))
  os.system("cat head origtxyz tmp >%s_EP.txyz"%(filename))

  '''remove the temp files'''
  os.system("rm -f tmp head origtxyz %s.xyz" %filename)
  del mol
  return 

main(sys.argv[1])

