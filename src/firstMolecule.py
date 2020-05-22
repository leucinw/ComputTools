
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# find the number of atoms of the first molecule in a molecular cluster

import os

def main(xyz):
  convertstr = "babel -ixyz %s -otxyz %s"  %(xyz, xyz.replace(".xyz", ".txyz"))
  os.system(convertstr)
  connections = []
  for line in open(xyz.replace(".xyz", ".txyz")).readlines()[1:]:
    d = line.split()
    connections.append([d[0]]+d[6:])

  firstmoleculeatoms = connections[0]
  connections[0] = []
  
  while True:
    precount = connections.count([])
    for i in range(1,len(connections)):
      for j in connections[i]: 
        if j in firstmoleculeatoms:
          firstmoleculeatoms += connections[i]
          connections[i] = []
    postcount = connections.count([])
    if precount == postcount:break
  return connections.count([])

if __name__ == "__main__":
  main()
