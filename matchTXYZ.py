

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Convert the second tinker xyz file to the first one

# $ python matchTXYZ.py template.xyz dealwith.xyz 
# result a dealwith.xyz_2 file

import os,sys

def readTXYZ(TXYZ, singleAtom=False, tailOnly=False):
  lines = open(TXYZ).readlines()[1:] 
  atoms=[];coord=[]
  order=[];types=[];connections=[]
  tails = []
  for line in lines:
    data=line.split()
    order.append(data[0])
    types.append(data[5])
    connections.append(data[6:])
    idxstr = " " + data[5] + " "
    idx = line.index(idxstr)
    tails.append(line[idx:])
    if singleAtom:
      atoms.append(data[1][0])
    else:
      atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
  if tailOnly:
    return tails
  else:
    return atoms,coord,order,types,connections

def fingerprint(TXYZ):
  atoms, _, order, _, connections =  readTXYZ(TXYZ, singleAtom=True)
  atom_dict = {}
  fp = []
  for i, j in zip(order, connections):
    if i not in atom_dict:
      atom_dict[i] = j 
  for i,j in atom_dict.items():
    mystr = ''
    for k in j:
      l = atom_dict[k]
      for m in l:
        n = atom_dict[m]
        for o in n:
          mystr += atoms[int(o)-1]
    mystr = ''.join(sorted(mystr)) 
    mystr = atoms[int(i)-1] + "->" + mystr
    fp.append(mystr)
  return fp


def main():
  template = sys.argv[1]
  dealwith = sys.argv[2]
  fname = dealwith + "_2" 
  fp1 = fingerprint(template)
  fp2 = fingerprint(dealwith)
  newidx = []
  for i in fp1:
    if i in fp2:
      idx = fp2.index(i)
      newidx.append(idx)
      fp2[idx] = ' '
    else:
      print("TWO STRUCTURES DIFFER! PLEASE CHECK !")
  atoms, coord, _, _, _ =  readTXYZ(dealwith, singleAtom=True)
  tails = readTXYZ(template,tailOnly=True)
  f = open(fname, "w")
  f.write("%3s\n"%len(atoms))
  for i in range(len(newidx)):
    idx = int(newidx[i])
    f.write("%3s%3s%12.6f%12.6f%12.6f   %s"%(i+1,atoms[idx], coord[idx][0], coord[idx][1], coord[idx][2], tails[i]))
    f.close()

if __name__ == "__main__":
  main()
