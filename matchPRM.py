

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Generate tinker prm file for new atomtype txyz

''' python matchPRM.py template.prm template.txyz dealwith.txyz energyterm'''

import os,sys,subprocess

tinkerexedir = "/home/liuchw/Softwares/tinkers/tinker8.2-intel13/source"

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
    if singleAtom:
      atoms.append(data[1][0])
    else:
      atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
  if tailOnly:
    return tails
  else:
    return atoms,coord,order,types,connections

def matchtorsion(templateprm, templatexyz, dealwithxyz):
  _,_,_,types1,_ = readTXYZ(templatexyz)
  _,_,_,types2,_ = readTXYZ(dealwithxyz)
  newidx = [i for i in range(len(types1))]
  typesDict = {}
  ttt = open("torsion", "w")
  for i in range(len(newidx)):
    if types1[i] not in typesDict:
      typesDict[types1[i]] = types2[int(newidx[i])]
  with open("tinker.key",'w') as f:
    f.write("parameters %s\n"%templateprm)
    f.write("torsionterm only\n")
  cmdstr = "%s/analyze.x %s -k tinker.key P > tmp.out "%(tinkerexedir, templatexyz)   
  subprocess.run(cmdstr,shell=True)
  lines = open("tmp.out").readlines()
  for line in lines:
    if "Torsional Angle Parameters :" in line:
      idx = lines.index(line)
  idx +=4
   
  torsionDict = {}
  onefold   = "0.000    0.0   1   "
  twofold   = "0.000  180.0   2   "
  threefold = "0.000    0.0   3   "
  repeated = []
  for i in range(idx, len(lines), 1):
    a1, a2, a3, a4 = lines[i].split()[1:5]
    if "0/1" in lines[i]:
      index = int(lines[i].split().index("0/1")) - 1
      onefold_ = '   '.join([lines[i].split()[index], '0.0', '1'])
    else:
      onefold_ = onefold
    if "180/2" in lines[i]:
      index = int(lines[i].split().index("180/2")) - 1
      twofold_ = '   '.join([lines[i].split()[index], '180.0', '2'])
    else:
      twofold_ = twofold
    if "0/3" in lines[i]:
      index = int(lines[i].split().index("0/3")) - 1
      threefold_ = '   '.join([lines[i].split()[index], '0.0', '3'])
    else:
      threefold_ = threefold
    torprm = ' '.join([onefold_, twofold_, threefold_])
    a1_ = typesDict[types1[int(a1)-1]]
    a2_ = typesDict[types1[int(a2)-1]]
    a3_ = typesDict[types1[int(a3)-1]]
    a4_ = typesDict[types1[int(a4)-1]]
    torstring = 4*"%5s"%(a1_, a2_, a3_, a4_)
    if torstring not in repeated:
      ttt.write("%s%5s%5s%5s%5s   %s\n"%("torsion", a1_, a2_, a3_, a4_, torprm))
      repeated.append(torstring)
  ttt.close()
  return

def matchbond(templateprm, templatexyz, dealwithxyz):
  _,_,_,types1,_ = readTXYZ(templatexyz)
  _,_,_,types2,_ = readTXYZ(dealwithxyz)
  newidx = [i for i in range(len(types1))]
  typesDict = {}
  ttt = open("bond", "w")
  for i in range(len(newidx)):
    if types1[i] not in typesDict:
      typesDict[types1[i]] = types2[int(newidx[i])]
  with open("tinker.key",'w') as f:
    f.write("parameters %s\n"%templateprm)
    f.write("bondterm only\n")
  cmdstr = "%s/analyze.x %s -k tinker.key P > tmp.out "%(tinkerexedir, templatexyz)   
  subprocess.run(cmdstr,shell=True)
  lines = open("tmp.out").readlines()
  for line in lines:
    if "Bond Stretching Parameters :" in line:
      idx = lines.index(line)
  idx +=4
   
  bondDict = {}
  repeated = []
  for i in range(idx, len(lines), 1):
    a1, a2, ks, bnd = lines[i].split()[1:5]
    bndprm = ' '.join([ks, bnd])
    a1_ = typesDict[types1[int(a1)-1]]
    a2_ = typesDict[types1[int(a2)-1]]
    bndstring = 2*"%5s"%(a1_, a2_)
    if bndstring not in repeated:
      ttt.write("%s%5s%5s  %s\n"%("bond", a1_, a2_, bndprm))
      repeated.append(bndstring)
  ttt.close()
  return

def matchangle(templateprm, templatexyz, dealwithxyz):
  _,_,_,types1,_ = readTXYZ(templatexyz)
  _,_,_,types2,_ = readTXYZ(dealwithxyz)
  newidx = [i for i in range(len(types1))]
  typesDict = {}
  ttt = open("angle", "w")
  for i in range(len(newidx)):
    if types1[i] not in typesDict:
      typesDict[types1[i]] = types2[int(newidx[i])]
  with open("tinker.key",'w') as f:
    f.write("parameters %s\n"%templateprm)
    f.write("angleterm only\n")
  cmdstr = "%s/analyze.x %s -k tinker.key P > tmp.out "%(tinkerexedir, templatexyz)   
  subprocess.run(cmdstr,shell=True)
  lines = open("tmp.out").readlines()
  for line in lines:
    if "Angle Bending Parameters :" in line:
      idx = lines.index(line)
  idx +=4
   
  angleDict = {}
  repeated = []
  for i in range(idx, len(lines), 1):
    a1, a2, a3, kb, ang = lines[i].split()[1:6]
    angprm = ' '.join([kb, ang])
    a1_ = typesDict[types1[int(a1)-1]]
    a2_ = typesDict[types1[int(a2)-1]]
    a3_ = typesDict[types1[int(a3)-1]]
    angstring = 3*"%5s"%(a1_, a2_, a3_)
    if angstring not in repeated:
      ttt.write("%s%5s%5s%5s  %s\n"%("angle", a1_, a2_, a3_, angprm))
      repeated.append(angstring)
  ttt.close()
  return

def matchopbend(templateprm, templatexyz, dealwithxyz):
  _,_,_,types1,_ = readTXYZ(templatexyz)
  _,_,_,types2,_ = readTXYZ(dealwithxyz)
  newidx = [i for i in range(len(types1))]
  typesDict = {}
  ttt = open("opbend", "w")
  for i in range(len(newidx)):
    if types1[i] not in typesDict:
      typesDict[types1[i]] = types2[int(newidx[i])]
  with open("tinker.key",'w') as f:
    f.write("parameters %s\n"%templateprm)
    f.write("opbendterm only\n")
  cmdstr = "%s/analyze.x %s -k tinker.key P > tmp.out "%(tinkerexedir, templatexyz)   
  subprocess.run(cmdstr,shell=True)
  lines = open("tmp.out").readlines()
  for line in lines:
    if "Out-of-Plane Bend Parameters :" in line:
      idx = lines.index(line)
  idx +=4
   
  opbendDict = {}
  repeated = []
  for i in range(idx, len(lines), 1):
    a1, a2, a3, a4, kopb = lines[i].split()[1:6]
    opbendprm = ' '.join([kopb])
    a1_ = typesDict[types1[int(a1)-1]]
    a2_ = typesDict[types1[int(a2)-1]]
    a3_ = typesDict[types1[int(a3)-1]]
    a4_ = typesDict[types1[int(a4)-1]]
    opbendstring = 4*"%5s"%(a1_, a2_, a3_, a4_)
    if opbendstring not in repeated:
      ttt.write("%s%5s%5s%5s%5s  %s\n"%("opbend", a1_, a2_, a3_, a4_, opbendprm))
      repeated.append(opbendstring)
  ttt.close()
  return
  
def main():
  templateprm = sys.argv[1]
  templatexyz = sys.argv[2]
  dealwithxyz = sys.argv[3]
  energyterm  = sys.argv[4]

  if energyterm.upper() == "TORSION":
    matchtorsion(templateprm, templatexyz, dealwithxyz)
  if energyterm.upper() == "BOND":
    matchbond(templateprm, templatexyz, dealwithxyz)
  if energyterm.upper() == "ANGLE":
    matchangle(templateprm, templatexyz, dealwithxyz)
  if energyterm.upper() == "OPBEND":
    matchopbend(templateprm, templatexyz, dealwithxyz)
     
  return

if __name__ == "__main__":
  main()
