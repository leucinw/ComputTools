
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

from pybel import *
import argparse
import numpy as np
import sys
# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

def genAtomType(txyz, potent):
  fname, _ = os.path.splitext(txyz)
  lines = open(txyz).readlines()
  if len(lines[0].split()) == 1:
    with open(txyz, "w") as f:
      f.write(lines[0].split("\n")[0] + " comments\n")
      for i in range(1,len(lines)):
        f.write(lines[i])
  types = np.loadtxt(txyz, usecols=(5,), unpack=True, dtype="str", skiprows=1)
  for mol in readfile('txyz',txyz):
    matchDict = {}
    matchList = []
    commentsDict = {} 
    classesDict = {} 
    natoms = len(mol.atoms)
    for i in range(1, natoms+1, 1):
      matchList.append(i)
    matchDict = dict.fromkeys(matchList, 0)
    commentsDict = dict.fromkeys(matchList, 0)
    classesDict = dict.fromkeys(matchList, 0)
    if potent.upper() == "POLAR":
      lines = open(os.path.join(typedir,"amoebaplusPolarType.dat")).readlines()
    elif potent.upper() == "CF":
      lines = open(os.path.join(typedir,"amoebaplusCFluxType.dat")).readlines()
    elif potent.upper() == "BONDED":
      lines = open(os.path.join(typedir,"amoebaplusBondedType.dat")).readlines()
    elif potent.upper() == "NONBONDED":
      lines = open(os.path.join(typedir,"amoebaplusNonbondedType.dat")).readlines()
    else:
      sys.exit(RED + f"{potent} not supported" + ENDC) 

    for line in lines:
      if "#" not in line[0] and len(line)>10:
        data = line.split()
        myStr = data[0]
        classNum = data[1]
        className = data[2]
        comment = line.split("# ")[1][0:-1]
        smarts = Smarts(myStr)
        match = smarts.findall(mol)
        if match:
          for i in range(len(match)):	
            matchDict[match[i][0]] = className
            commentsDict[match[i][0]] = comment
            classesDict[match[i][0]] = classNum
    with open(f"{fname}.type", "w") as f:	
      for i in range(1, natoms+1, 1):
        f.write("%5s %5s %5s %5s #%s\n"%(i, types[i-1], matchDict[i], classesDict[i], commentsDict[i]))
  return

def assignPolar(fname, tinkerkey):
  types, polars = np.loadtxt(os.path.join(prmdir,"polarize.prm"), usecols=(0,1), unpack=True, dtype="str",skiprows=1)
  smartspolarDict = dict(zip(types, polars))
  ttypes, stypes = np.loadtxt(f"{fname}.type", usecols=(1,2), unpack=True, dtype="str")
  tinkerpolarDict = {}
  for t,s in zip(ttypes, stypes): 
    if t not in tinkerpolarDict:
      tinkerpolarDict[t] = smartspolarDict[s]
  lines = open(tinkerkey).readlines()
  with open(tinkerkey, "w") as f:
    for line in lines:
      if "polarize " in line:
        dd = line.split()
        if dd[1] in ttypes:
          dd[2] = tinkerpolarDict[dd[1]]
          line = "    ".join(dd) + "\n"
          print(GREEN + "polarizability parameter found for %s"%dd[1] + ENDC)
      f.write(line)
  return

def assignNonbonded(fname, tinkerkey):
  types = np.loadtxt(os.path.join(prmdir,"nonbonded.prm"), usecols=(0,), unpack=True, dtype="str",skiprows=2)
  cpalphas, cpnucs = np.loadtxt(os.path.join(prmdir,"nonbonded.prm"), usecols=(1,2), unpack=True, dtype="float",skiprows=2)
  ctas, ctbs = np.loadtxt(os.path.join(prmdir,"nonbonded.prm"), usecols=(3,4), unpack=True, dtype="float",skiprows=2)
  vdwrs, vdwes, vdwreds = np.loadtxt(os.path.join(prmdir,"nonbonded.prm"), usecols=(5,6,7), unpack=True, dtype="float",skiprows=2)
  CP = [[i,j] for i,j in zip(cpalphas, cpnucs)]
  CT = [[i,j] for i,j in zip(ctas, ctbs)]
  VDW = [[i,j,k] for i,j,k in zip(vdwrs, vdwes, vdwreds)]
  smartsCPDict = dict(zip(types, CP))
  smartsCTDict = dict(zip(types, CT))
  smartsVDWDict = dict(zip(types, VDW))
  ttypes, stypes = np.loadtxt(f"{fname}.type", usecols=(1,2), unpack=True, dtype="str")
  tinkerCPDict = {}
  tinkerCTDict = {}
  tinkerVDWDict = {}
  for t,s in zip(ttypes, stypes): 
    if t not in tinkerCPDict:
      tinkerCPDict[t] = smartsCPDict[s]
      tinkerCTDict[t] = smartsCTDict[s]
      tinkerVDWDict[t] = smartsVDWDict[s]
  lines = open(tinkerkey).readlines()
  with open(tinkerkey, "a") as f:
    f.write("# charge penetration parameters assigned from database\n")
    for t in tinkerCPDict:
      ''' the commented out line is for old tinker 8.2 '''
      #line = "cp  %5s%10.5f%10.5f\n"%(t, tinkerCPDict[t][0], tinkerCPDict[t][1])
      line = "chgpen  %5s%10.5f%10.5f\n"%(t, tinkerCPDict[t][1], tinkerCPDict[t][0])
      f.write(line)
    print(GREEN + "charge penetration parameters assigned from database"+ ENDC)
    f.write("# charge transfer parameters assigned from database\n")
    for t in tinkerCTDict:
      #line = "ct  %5s%10.5f%10.5f\n"%(t, tinkerCTDict[t][0], tinkerCTDict[t][1])
      line = "chgtrn  %5s%10.5f%10.5f\n"%(t, tinkerCTDict[t][0]*3.01147, tinkerCTDict[t][1])
      f.write(line)
    print(GREEN + "charge transfer parameters assigned from database" + ENDC)
    f.write("# van der Waals parameters assigned from database\n")
    for t in tinkerVDWDict:
      if tinkerVDWDict[t][2] != 1.0:
        line = "vdw  %5s%10.5f%10.5f%10.5f\n"%(t, tinkerVDWDict[t][0], tinkerVDWDict[t][1], tinkerVDWDict[t][2])
      else:
        line = "vdw  %5s%10.5f%10.5f\n"%(t, tinkerVDWDict[t][0], tinkerVDWDict[t][1])
      f.write(line)
    print(GREEN+"van der Waals parameters assigned from database"+ENDC)
  return

def assignCFlux(fname, tinkerkey):
  # map from tinker type number to smart type string
  # this will be used by bndcflux and angcflux
  ttypes, stypes = np.loadtxt(f"{fname}.type", usecols=(1,2), unpack=True, dtype="str")
  tinker2smarts = dict(zip(ttypes, stypes))
  #cflux-b
  '''bond cflux atom indices are interchangable'''
  type1, type2, jbonds = np.loadtxt(os.path.join(prmdir,"cfbond.prm"), usecols=(-2,-1,1), unpack=True, dtype="str",skiprows=1)
  types = []
  for t1, t2 in zip(type1, type2):
    types.append(t1 + "_" + t2)
  smartsCFbondDict = dict(zip(types, jbonds))
  lines = open(tinkerkey).readlines()

  with open(tinkerkey,"a") as f:
    f.write("# CHGFLX parameters assigned from database\n")
    for line in lines:
      if "bond " in line:
        dd = line.split()
        if (dd[1] in ttypes) and (dd[2] in ttypes):
          s1 = tinker2smarts[dd[1]]
          s2 = tinker2smarts[dd[2]]
          comb1 = s1 + "_" + s2
          comb2 = s2 + "_" + s1
          if comb1 in smartsCFbondDict:
            f.write("bndcflux %s %s %10.5f\n"%(dd[1], dd[2], float(smartsCFbondDict[comb1])))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          elif comb2 in smartsCFbondDict:
            f.write("bndcflux %s %s %10.5f\n"%(dd[1], dd[2], float(smartsCFbondDict[comb2])))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          else:
            print(RED + "CFlux parameter NOT found for bond %s-%s"%(dd[1], dd[2]) + ENDC)

  #cflux-a
  '''angle cflux in parameter file is in the right order for jt1,jt2,jb1,jb2'''
  '''when assign parameters, need to first sort the angle atom indices, then to match database'''
  type1, type2, type3  = np.loadtxt(os.path.join(prmdir,"cfangle.prm"), usecols=(-3, -2, -1), unpack=True, dtype="str",skiprows=1)
  jtheta1, jtheta2, jbond1, jbond2  = np.loadtxt(os.path.join(prmdir,"cfangle.prm"), usecols=(1, 2, 3, 4), unpack=True, dtype="float",skiprows=1)

  types = []
  jparams = []

  '''store two sets of parameters considering the assymetry of angle-cflux'''
  for t1, t2, t3 in zip(type1, type2, type3):
    types.append(t1 + "_" + t2 + "_" + t3)
    types.append(t3 + "_" + t2 + "_" + t1)

  for jt1, jt2, jb1, jb2 in zip(jtheta1, jtheta2, jbond1, jbond2):
    # convert jt unit from e/degree to e/radian
    '''this will be compatitable with new tinker version >= 8.7'''
    jt1 *= 57.2958
    jt2 *= 57.2958
    jparams.append(" ".join(["%10.5f"%jt1, "%10.5f"%jt2, "%10.5f"%jb1, "%10.5f"%jb2]))
    jparams.append(" ".join(["%10.5f"%jt2, "%10.5f"%jt1, "%10.5f"%jb2, "%10.5f"%jb1]))
  
  smartsCFangleDict = dict(zip(types, jparams))
  
  with open(tinkerkey, "a") as f:
    for line in lines:
      if "angle " in line:
        dd = line.split()
        angletype1 = dd[1]
        angletype2 = dd[2]
        angletype3 = dd[3]
        if (angletype1 in ttypes) and (angletype2 in ttypes) and (angletype3 in ttypes):
          '''always make sure type1 <= type3'''
          if int(angletype1) > int(angletype3):
            angletype1, angletype3 = angletype3, angletype1
            print("flipped angletype1 and angletype3")
          s1 = tinker2smarts[angletype1]
          s2 = tinker2smarts[angletype2]
          s3 = tinker2smarts[angletype3]
          comb = s1 + "_" + s2 + "_" + s3
          if comb in smartsCFangleDict:
            f.write("angcflux %s %s %s %s\n"%(angletype1, angletype2, angletype3, smartsCFangleDict[comb]))
            print(GREEN + "CFlux parameters found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          else:
            print(RED + "CFlux parameters NOT found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
  return

def assignBonded(fname, tinkerkey):
  # map from atom type to atom class
  # this will be used by all the following bonded terms
  ttypes, tclasses = np.loadtxt(f"{fname}.type", usecols=(1,3), unpack=True, dtype="str")
  type2class = dict(zip(ttypes, tclasses))

  # bond stretching
  # only assign force constant parameter, since equilibrium length will be from QM
  class1, class2, kbonds = np.loadtxt(os.path.join(prmdir,"bond.prm"), usecols=(1,2,3), unpack=True, dtype="str",skiprows=1)
  classes = []
  for c1, c2 in zip(class1, class2):
    classes.append(c1 + "_" + c2)
  classKconstantDict = dict(zip(classes, kbonds))
  lines = open(tinkerkey).readlines()
  with open(tinkerkey,"a") as f:
    f.write("# BONDED parameters assigned from database\n")
    for line in lines:
      if "bond " in line:
        dd = line.split()
        if (dd[1] in ttypes) and (dd[2] in ttypes):
          c1 = type2class[dd[1]]
          c2 = type2class[dd[2]]
          comb1 = c1 + "_" + c2
          comb2 = c2 + "_" + c1
          if comb1 in classKconstantDict:
            f.write("bond %s %s %10.5f %s\n"%(dd[1], dd[2], float(classKconstantDict[comb1]), dd[4]))
            print(GREEN + "BOND stretching parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          elif comb2 in classKconstantDict:
            f.write("bond %s %s %10.5f %s\n"%(dd[1], dd[2], float(classKconstantDict[comb2]), dd[4]))
            print(GREEN + "BOND stretching parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          else:
            print(RED + "BOND stretching parameter NOT found for bond %s-%s with atom class %s"%(dd[1], dd[2], comb1) + ENDC)

  #angle bending
  class1, class2, class3  = np.loadtxt(os.path.join(prmdir, "angle.prm"), usecols=(1, 2, 3), unpack=True, dtype="str",skiprows=1)
  angleKconstant  = np.loadtxt(os.path.join(prmdir, "angle.prm"), usecols=(4), unpack=True, dtype="float",skiprows=1)
  classes = []
  angKconsts = []
  '''store two sets of parameters since angle indices are interchangable''' 
  for c1, c2, c3 in zip(class1, class2, class3):
    classes.append(c1 + "_" + c2 + "_" + c3)
    classes.append(c3 + "_" + c2 + "_" + c1)
  for k in angleKconstant:
    angKconsts.append(k)
    angKconsts.append(k)

  classAngleKconstantDict = dict(zip(classes, angKconsts))
  
  with open(tinkerkey, "a") as f:
    for line in lines:
      if "angle " in line:
        dd = line.split()
        angletype1 = dd[1]
        angletype2 = dd[2]
        angletype3 = dd[3]
        if (angletype1 in ttypes) and (angletype2 in ttypes) and (angletype3 in ttypes):
          c1 = type2class[angletype1]
          c2 = type2class[angletype2]
          c3 = type2class[angletype3]
          comb = c1 + "_" + c2 + "_" + c3
          if (comb in classAngleKconstantDict): 
            f.write("angle %s %s %s %s %s\n"%(angletype1, angletype2, angletype3, classAngleKconstantDict[comb], dd[5]))
            print(GREEN + "ANGLE bending parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          else: 
            print(RED + "ANGLE bending parameter NOT found for angle %s-%s-%s with atom class %s"%(angletype1, angletype2, angletype3, comb) + ENDC)

  #bond-angle coupling (strbnd term)
  class1, class2, class3  = np.loadtxt(os.path.join(prmdir, "strbend.prm"), usecols=(1, 2, 3), unpack=True, dtype="str",skiprows=1)
  strbndK1, strbndK2  = np.loadtxt(os.path.join(prmdir, "strbend.prm"), usecols=(4,5), unpack=True, dtype="float",skiprows=1)
  classes = []
  strbndKs = []
  '''store two sets of parameters since strbnd is asymetric''' 
  for c1, c2, c3 in zip(class1, class2, class3):
    classes.append(c1 + "_" + c2 + "_" + c3)
    classes.append(c3 + "_" + c2 + "_" + c1)
  for k1, k2 in zip(strbndK1, strbndK2):
    strbndKs.append("%10.4f%10.4f"%(k1,k2))
    strbndKs.append("%10.4f%10.4f"%(k2,k1))
  classStrbndKconstantDict = dict(zip(classes, strbndKs))
  
  with open(tinkerkey, "a") as f:
    for line in lines:
      #if "strbnd " in line:
      # try to find strbnd parameters for every angle
      if "angle " in line:
        dd = line.split()
        angletype1 = dd[1]
        angletype2 = dd[2]
        angletype3 = dd[3]
        if (angletype1 in ttypes) and (angletype2 in ttypes) and (angletype3 in ttypes):
          c1 = type2class[angletype1]
          c2 = type2class[angletype2]
          c3 = type2class[angletype3]
          comb = c1 + "_" + c2 + "_" + c3
          if (comb in classStrbndKconstantDict):
            f.write("strbnd %s %s %s %s\n"%(angletype1, angletype2, angletype3, classStrbndKconstantDict[comb]))
            print(GREEN + "STRBND coupling parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          else: 
            print(RED + "STRBND coupling parameter NOT found for angle %s-%s-%s with atom class %s"%(angletype1, angletype2, angletype3, comb) + ENDC)

  # out-of-plane bending 
  class1, class2, kbonds = np.loadtxt(os.path.join(prmdir,"opbend.prm"), usecols=(1,2,5), unpack=True, dtype="str",skiprows=1)
  classes = []
  for c1, c2 in zip(class1, class2):
    classes.append(c1 + "_" + c2)
  classKconstantDict = dict(zip(classes, kbonds))
  lines = open(tinkerkey).readlines()
  with open(tinkerkey,"a") as f:
    for line in lines:
      if "opbend " in line:
        dd = line.split()
        if (dd[1] in ttypes) and (dd[2] in ttypes):
          c1 = type2class[dd[1]]
          c2 = type2class[dd[2]]
          comb = c1 + "_" + c2
          if comb1 in classKconstantDict:
            f.write("opbend %s %s    0    0 %s\n"%(dd[1], dd[2], float(classKconstantDict[comb])))
            print(GREEN + "OPBEND parameter assigned for bond %s-%s-0-0"%(dd[1], dd[2]) + ENDC)
          else:
            print(RED + "OPBEND parameter NOT found for bond %s-%s-0-0 with atom class %s"%(dd[1], dd[2], comb) + ENDC)
  return
 
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'xyz', help = "tinker xyz file", required=True)  
  parser.add_argument('-k', dest = 'key', help = "tinker prm file", required=True)  
  parser.add_argument('-p', dest = 'potent', help = "potential energy term", required=True, choices=["POLAR", "CF", "BONDED", "NONBONDED"])  
  args = vars(parser.parse_args())
  inp = args["xyz"]
  key = args["key"]
  potent = args["potent"]
  global typedir 
  global prmdir 
  typedir="/home/liuchw/amoebaplus.data/types"
  prmdir="/home/liuchw/amoebaplus.data/prms"
  genAtomType(inp, potent)
  fname, _ = os.path.splitext(inp)
  if (potent.upper() == "POLAR"):
    assignPolar(fname, key)
  elif (potent.upper() == "CF"):
    assignCFlux(fname, key) 
  elif (potent.upper() == "BONDED"):
    assignBonded(fname, key) 
  elif (potent.upper() == "NONBONDED"):
    assignNonbonded(fname, key) 
  else:
    print(RED + f"{potent} term not supported!" + ENDC)
  return

if len(sys.argv) == 1:
  print('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
else:
  main()
