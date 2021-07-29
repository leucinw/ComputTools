
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


''' Usage: python matchTXYZ.py template.txyz dealwith.(t)xyz 
    # Assign the atom types of the template.txyz file to the dealwith txyz or xyz file.
    # Especially suitable for the case that their atoms are different in order.
'''

import os
import sys
import numpy as np

def readTXYZ(TXYZ):
  atoms=[]
  coord=[]
  order=[]
  types=[]
  connections=[]
  for line in open(TXYZ).readlines()[1:]: 
    data=line.split()
    order.append(data[0])
    atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
    types.append(data[5])
    connections.append(data[6:])
  return atoms,coord,order,types,connections

def fingerprint(TXYZ):
  fprints = []
  atoms, elements = np.loadtxt(TXYZ, usecols=(0,1), dtype='str', skiprows=1, unpack=True)
  connections = []
  for line in open(TXYZ).readlines()[1:]:
    d = line.split()
    connections.append(d[6:])
  
  atom_ele_dict = dict(zip(atoms, elements))
  atom_con_dict = {}
  for atom, con in zip(atoms,connections):
    con_ele = [atom_ele_dict[c] for c in con] 
    constr = ''.join(sorted(con_ele)) 
    atom_con_dict[atom] = constr

  level = 5 
  if level > 1:
    atom_con_dict2 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict[c] for c in cons])
      atom_con_dict2[atom] = ''.join(sorted(newstr))

  # level 3 is good for chain molecules 
  if level > 2:
    atom_con_dict3 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict2[c] for c in cons])
      atom_con_dict3[atom] = ''.join(sorted(newstr))

  # level 4 is needed for ring molecules 
  if level > 3:
    atom_con_dict4 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict3[c] for c in cons])
      atom_con_dict4[atom] = ''.join(sorted(newstr))
  
  if level > 4:
    atom_con_dict5 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict4[c] for c in cons])
      atom_con_dict5[atom] = ''.join(sorted(newstr))
  
  for atom in atoms:
    fprints.append(atom_ele_dict[atom] + '-' + str(''.join(sorted(atom_con_dict[atom] + atom_con_dict2[atom] + atom_con_dict3[atom] + atom_con_dict4[atom] + atom_con_dict5[atom]))))
  
  return fprints

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
  for i in fp2:
    if i in fp1:
      idx = fp1.index(i)
      newidx.append(idx)
      fp1[idx] = ' '
    else:
      print(f"Error:{template} and {dealwith} could not match!")
  
  atoms, coord, _, _, connections =  readTXYZ(dealwith)
  _, _, _,types, _ =  readTXYZ(template)
  with open(fname, 'w') as f:
    f.write("%3s\n"%len(atoms))
    for i in range(len(newidx)):
      idx = int(newidx[i])
      f.write("%3s%3s%12.6f%12.6f%12.6f  %s   %s\n"%(i+1,atoms[i], coord[i][0], coord[i][1], coord[i][2], types[idx], '  '.join(connections[i])))
  return 

if __name__ == "__main__":
  main()
