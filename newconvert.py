
from readChemFile import *
from writeChemFile import *
import os

def ToXYZ(f):
  fname, ext = os.path.splitext(f)
  fname += ".xyz"
  if ext == ".log":
    atoms, coords, _ = readLOG(f)
  elif ext == ".txyz":
    atoms, coords, _,_,_ = readTXYZ(f)
  elif ext == ".com":
    atoms, coords, _ = readCOM(f)
  elif ext == ".xyz":
    atoms, coords= readXYZ(f)
    fname += "_2"
  else:
    print("Error: the input file type not supported")
  if atoms != []:
    writeXYZ(fname, atoms, coords)
  return

def ToTXYZ(f):
  fname, ext = os.path.splitext(f)
  def getConnection(xyz):
    orders = []
    connections = []
    os.system("babel -ixyz %s -opdb tmp.pdb"%xyz)
    for line in open("tmp.pdb").readlines():
      if "CONECT" in line:
        d = line.split()
        orders.append(d[1])
        connections.append(d[2:])
    return orders,connections
  if ext == ".txyz":
    fname = f + "_2"
    atoms, coords, orders, types, connections = readTXYZ(f)
    writeTXYZ(fname, orders, atoms, coords, types, connections)
  else:
    ToXYZ(f)
    xyz = fname + ".xyz"
    txyz = fname + ".txyz"
    os.system("babel -ixyz %s -otxyz %s"%(xyz, txyz))
  return


import sys
ToTXYZ(sys.argv[1])
