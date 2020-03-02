

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Using fingerprint to automatically match the atom type, this script has two functions:
# 1. Generate a tinker-xyz file based on a "normal xyz" file and a template tinker-xyz file
# 2. Generate a tinker-xyz file based on a "tinker xyz" file and a template tinker-xyz file

''' python matchTXYZ.py template.txyz dealwith.(t)xyz '''

import os,sys
#========================#
# read tinker-xyz file   #
#========================#
def readTXYZ(TXYZ, singleAtom=False, tailOnly=False):
  atoms=[];coord=[]; tails = []
  order=[];types=[];connections=[]
  for line in open(TXYZ).readlines()[1:]: 
    data=line.split()
    order.append(data[0])
    types.append(data[5])
    connections.append(data[6:])
    idxstr = " " + data[5] + " "
    idx = line.index(idxstr)
    tails.append(line[idx:])
    #=================================================
    # sometime AA is "CA", it is Carbon, not Calcium 
    # ================================================
    if singleAtom:
      atoms.append(data[1][0])
    else:
      atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
  #=======================================#
  # return only connection and type info. #
  #=======================================#
  if tailOnly:
    return tails
  #=======================================#
  # return all info. of a tinker-xyz file #
  #=======================================#
  else:
    return atoms,coord,order,types,connections

#=============================================
# use .mna provided by obabel as fingerprint
#=============================================
def fingerprint(TXYZ):
  #===========================================#
  # add some comments if no comments in txyz  #
  #===========================================#
  lines = open(TXYZ).readlines()
  if len(lines[0].split()) == 1:
    with open(TXYZ, 'w') as f:
      f.write(lines[0][:-1] + " Comments required by obabel\n")
      for line in lines[1:]:
        f.write(line)
  fp = []
  MNA = TXYZ.replace("txyz", "mna")
  obstr = "obabel -itxyz %s -omna -O %s"%(TXYZ, MNA)
  print(obstr)
  os.system(obstr)
  for line in open(MNA).readlines()[2:]:
    fpstr = ''.join(sorted(line[:-1])) 
    fp.append(fpstr)
  return fp

#=======================================================#
# make a xyz2txyz convertion if normal xyz is provided  #
#=======================================================#
def xyz2txyz(xyz):
  txyz = xyz.replace("xyz", "txyz")
  obstr = "obabel -ixyz %s -otxyz -O %s"%(xyz, txyz)
  os.system(obstr)
  return

#===============#
# main function #
#===============#
def main():
  template = sys.argv[1]
  dealwith = sys.argv[2]
  if os.path.splitext(dealwith)[1] == ".xyz":
    xyz2txyz(dealwith)
    dealwith = dealwith.replace("xyz", "txyz")
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
  with open(fname, 'w') as f:
    f.write("%3s\n"%len(atoms))
    for i in range(len(newidx)):
      idx = int(newidx[i])
      f.write("%3s%3s%12.6f%12.6f%12.6f   %s"%(i+1,atoms[idx], coord[idx][0], coord[idx][1], coord[idx][2], tails[i]))
  return

if __name__ == "__main__":
  main()
