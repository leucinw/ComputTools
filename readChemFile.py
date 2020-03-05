
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


from readDataFile import *
from chem import *
# read tinker xyz file
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

# read gaussian log file
def readLOG(LOG):
  atoms=[]; data =[]; coord=[]; chargeMultip=[]
  with open(LOG) as f: 
    lines=f.readlines()
    for n in range(len(lines)):
      if "Charge =" in lines[n]:
        charge = (lines[n].split()[2])   
      if "Multiplicity =" in lines[n]:
        multip = (lines[n].split()[5])   
      if "NAtoms=" in lines[n]:
        natoms=int(lines[n].split()[1])   
      if ("Standard orientation" in lines[n]) or ("Input orientation" in lines[n]): 
        for m in range(5+n,5+natoms+n,1):
          data.append(lines[m][0:-1])
    chargeMultip = [charge, multip]  
    for i in range(natoms):
      d = data[i-natoms].split()
      atoms.append(pTable[d[1]])
      coord.append([float(d[3]), float(d[4]), float(d[5])])
  return atoms,coord,chargeMultip

# read gaussian input file
def readCOM(COM):
  atoms=[]; coord=[]; chargeMultip = []
  cpcorr = False
  for line in open(COM).readlines():
    if line[0] == "#":
      routeSec = line[:-1]
      upperRouteSec = routeSec.upper()
      if "COUNTERPOISE" in upperRouteSec:
        cpcorr = True
    if cpcorr == True:
      if "Fragment" in line:
        d = line.split()
        atoms.append(d[0].split("(")[0])
        coord.append([float(d[1]), float(d[2]), float(d[3])])
    else:
      if len(line.split()) == 4 and len(line) > 38 and (line.split()[0] in eTable):
        d = line.split()
        atoms.append(d[0])
        coord.append([float(d[1]), float(d[2]), float(d[3])])
  return atoms,coord,routeSec

# read common xyz file
def readXYZ(XYZ):
  atoms=[]; coord=[]
  for line in open(XYZ).readlines()[2:]:
    d = line.split()
    atoms.append(d[0])
    coord.append([float(d[1]), float(d[2]), float(d[3])])
  return atoms,coord

def readPsi4Out(OUT):
  lines = readWholeFile(OUT)
  ele  = 0.0
  exch = 0.0
  ind  = 0.0
  disp = 0.0
  tot  = 0.0
  for line in lines:
    if "Electrostatics     " in line:
      ele = float(line.split()[3])
    elif "Exchange     " in line:
      exch = float(line.split()[3])
    elif "Induction     " in line:
      ind = float(line.split()[3])
    elif "Dispersion     " in line:
      disp = float(line.split()[3])
  vdw = exch + disp
  tot = vdw + ind + ele
  return ele, ind, exch, disp, vdw, tot
  
def readTinkerOut(OUT):
  lines = readWholeFile(OUT)
  ele  = 0.0
  vdw = 0.0
  pol = 0.0
  tot  = 0.0
  for line in lines:
    if "Electrostatics     " in line:
      ele = float(line.split()[3])
    elif "Exchange     " in line:
      exch = float(line.split()[3])
    elif "Induction     " in line:
      ind = float(line.split()[3])
    elif "Dispersion     " in line:
      disp = float(line.split()[3])
  vdw = exch + disp
  tot = vdw + ind + ele
  return ele, ind, exch, disp, vdw, tot

# read tinker parameter file

def readTinkerPRM(PRM, term):
  lines = readWholeFile(PRM)
  # atom
  if term == "atom ":
    atoms = []
    atomList = []
    typeList = []
    remainder = []
    for line in lines:
      if term in line:
        data = line.split()
        atoms.append(data[0])
        atomList.append(int(data[1]))
        typeList.append(int(data[2]))
        remainder.append(data[3:])
    return atoms, atomList, typeList, remainder

  if term == "multipole ":
    mpoleList = []  
  if term == "polarize ":
    polarizes = []  
    atomList = []
    atomicPolarList = []
    dampingList = []
    polarGroupList = []
    for line in lines:
      if term in line:
        data = line.split()
        polarizes.append(data[0])
        atomList.append(int(data[1]))
        atomicPolarList.append(data[2])
        # good for AMOEBA not for AMOEBA+
        dampingList.append(data[3])
        polarGroupList.append([int(data[i]) for i in range(4, len(data), 1)])
    return polarizes, atomList, atomicPolarList, dampingList, polarGroupList

  if term == "vdw ":
    vdwList = []  
  if term == "bond ":
    bondList = []  
  if term == "angle ":
    angleList = []  
  if term == "strbnd ":
    strbndList = []
  return
#import sys
#print(readLOG(sys.argv[1]) )
#print(readCOM(sys.argv[1]) )
