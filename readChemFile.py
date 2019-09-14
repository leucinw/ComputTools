#!/usr/bin/env python

from readDataFile import *

# read tinker xyz file
def readTXYZ(TXYZ):
  lines = readWholeFile(TXYZ)[1:] #TINKER coordinate starts from second line
  atoms=[];coord=[]
  order=[];types=[];connections=[]
  for line in lines:
    data=line.split()
    order.append(data[0])
    types.append(data[5])
    connections.append(data[6:])
    atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
  return atoms,coord,order, types, connections

# read gaussian output file
def readLOG(LOG):
  data=[];coord=[]
  natoms=0
  reach_coordinate=0
  f=open(LOG,'r')
  while True:
    lines=f.readlines()
    if not lines:break
    for n in range(len(lines)):
      if "Charge =" in lines[n]:
        charge=(lines[n].split()[2])   
      if "Multiplicity =" in lines[n]:
        multip=(lines[n].split()[5])   
      if "NAtoms=" in lines[n]:
        natoms=int(lines[n].split()[1])   
      if "Standard orientation" in lines[n] or "Input orientation" in lines[n]: # Coordinates from "Standard orientation" part
        for m in range(5,5+natoms,1):
          data.append(lines[n+m][0:-1])
  for i in range(natoms):
     coord.append(data[i-natoms])
  f.close()
  return coord,charge,multip

# read gaussian input file
def readCOM(COM):
  eTable = ['H','He','Li','Be','B','C','N','O',\
                 'F','Ne','Na','Mg','Al','Si','P',\
                 'S','Cl','Ar','K','Ca','Ti','Fe',\
                 'Co','Ni','Cu','Zn',"Br", "I", "Cs"]
  atoms=[];coord=[]
  lines = readWholeFile(COM)
  for line in lines:
    if len(line.split())==4 and len(line)>50 and (line.split()[0] in eTable):
      data = line.split()
      atoms.append(data[0])
      coord.append([float(data[1]), float(data[2]), float(data[3])])
  return atoms,coord

# read common xyz file
def readXYZ(XYZ):
  atoms = []
  coord = []
  lines = open(XYZ).readlines()
  for line in lines[2:]:
    data = line.split()
    if len(data)==4:
      atoms.append(data[0])
      coord.append([float(data[1]), float(data[2]), float(data[3])])
  return atoms,coord

