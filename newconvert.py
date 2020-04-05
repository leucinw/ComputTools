import os


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
