

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# 1. Generate tinker-xyz by matching "normal-xyz to tinker-xyz"
# 2. Generate tinker-xyz by matching "tinker-xyz to tinker-xyz" 

''' python matchTXYZ.py template.txyz dealwith.(t)xyz '''

import os,sys

def readTXYZ(TXYZ, singleAtom=False, tailOnly=False):
  atoms=[];coord=[]; tails = []
  order=[];types=[];connections=[]
  for line in open(TXYZ).readlines()[1:]: 
    data=line.split()
    order.append(data[0])
    types.append(data[5])
    connections.append(data[6:])
    tails.append(' '.join(data[5:]))
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
  lines = open(TXYZ).readlines()
  if len(lines[0].split()) == 1:
    with open(TXYZ, 'w') as f:
      f.write(lines[0][:-1] + " Comments required by obabel\n")
      for line in lines[1:]:
        f.write(line)
  fp = []
  MNA = TXYZ.replace("txyz", "mna")
  obstr = "obabel -itxyz %s -omna -O %s"%(TXYZ, MNA)
  os.system(obstr)
  for line in open(MNA).readlines()[2:]:
    idx = line.index("(")
    atm = line[0:idx]
    fpstr = ''.join(sorted(line[idx:-1])) 
    fp.append(atm + fpstr)
  return fp

def main():
  template = sys.argv[1]
  dealwith = sys.argv[2]
  if os.path.splitext(dealwith)[1] == ".xyz":
    xyz = dealwith
    dealwith = dealwith.replace("xyz", "txyz")
    obstr = "obabel -ixyz %s -otxyz -O %s"%(xyz, dealwith)
    os.system(obstr)
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
      print("TWO STRUCTURES MAY BE DIFFERENT?")
  atoms, coord, _, _, _ =  readTXYZ(dealwith, singleAtom=True)
  tails = readTXYZ(template,tailOnly=True)
  with open(fname, 'w') as f:
    f.write("%3s\n"%len(atoms))
    for i in range(len(newidx)):
      idx = int(newidx[i])
      f.write("%3s%3s%12.6f%12.6f%12.6f   %s\n"%(i+1,atoms[idx], coord[idx][0], coord[idx][1], coord[idx][2], tails[i]))
  return

if __name__ == "__main__":
  main()
