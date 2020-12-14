
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
    natoms = len(mol.atoms)
    for i in range(1, natoms+1, 1):
      matchList.append(i)
    matchDict = dict.fromkeys(matchList, 0)
    commentsDict = dict.fromkeys(matchList, 0)
    if potent.upper() == "CF":
      lines = open(dbdir+"amoebaplusCFluxType.dat").readlines()
    elif potent.upper() == "POLAR":
      lines = open(dbdir+"amoebaplusPolarType.dat").readlines()
    elif potent.upper() == "NONBONDED":
      lines = open(dbdir+"amoebaplusNonbondedType.dat").readlines()
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
    with open(f"{fname}.type", "w") as f:	
      for i in range(1, natoms+1, 1):
        f.write("%3s %3s %3s #%s\n"%(i, types[i-1], matchDict[i], commentsDict[i]))
  return

def assignPolar(fname, tinkerkey):
  types, polars = np.loadtxt(dbdir+"polarize.prm", usecols=(0,1), unpack=True, dtype="str",skiprows=1)
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
  types = np.loadtxt(dbdir+"nonbonded.prm", usecols=(0,), unpack=True, dtype="str",skiprows=2)
  cpalphas, cpnucs = np.loadtxt(dbdir+"nonbonded.prm", usecols=(1,2), unpack=True, dtype="float",skiprows=2)
  ctas, ctbs = np.loadtxt(dbdir+"nonbonded.prm", usecols=(3,4), unpack=True, dtype="float",skiprows=2)
  vdwrs, vdwes, vdwreds = np.loadtxt(dbdir+"nonbonded.prm", usecols=(5,6,7), unpack=True, dtype="float",skiprows=2)
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
      line = "cp  %5s%10.5f%10.5f\n"%(t, tinkerCPDict[t][0], tinkerCPDict[t][1])
      f.write(line)
    print(GREEN + "charge penetration parameters assigned from database"+ ENDC)
    f.write("# charge transfer parameters assigned from database\n")
    for t in tinkerCTDict:
      line = "ct  %5s%10.5f%10.5f\n"%(t, tinkerCTDict[t][0], tinkerCTDict[t][1])
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
  #cflux-b
  type1, type2, jbonds = np.loadtxt(dbdir+"cfbond.prm", usecols=(-2,-1,1), unpack=True, dtype="str",skiprows=1)
  types = []
  for t1, t2 in zip(type1, type2):
    types.append(t1 + "_" + t2)
  smartsCFbondDict = dict(zip(types, jbonds))
  ttypes, stypes = np.loadtxt(f"{fname}.type", usecols=(1,2), unpack=True, dtype="str")
  tinker2smarts = dict(zip(ttypes, stypes))
  lines = open(tinkerkey).readlines()

  with open(tinkerkey,"a") as f:
    f.write("# cflux-b parameter assigned from database\n")
    for line in lines:
      if "bond " in line:
        dd = line.split()
        if (dd[1] in ttypes) and (dd[2] in ttypes):
          s1 = tinker2smarts[dd[1]]
          s2 = tinker2smarts[dd[2]]
          comb1 = s1 + "_" + s2
          comb2 = s2 + "_" + s1
          if comb1 in smartsCFbondDict:
            f.write("cflux-b %s %s %10.5f\n"%(dd[1], dd[2], float(smartsCFbondDict[comb1])))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          elif comb2 in smartsCFbondDict:
            f.write("cflux-b %s %s %10.5f\n"%(dd[1], dd[2], float(smartsCFbondDict[comb2])))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(dd[1], dd[2]) + ENDC)
          else:
            print(RED + "CFlux parameter NOT found for bond %s-%s"%(dd[1], dd[2]) + ENDC)

  #cflux-a
  type1, type2, type3  = np.loadtxt(dbdir+"cfangle.prm", usecols=(-3, -2, -1), unpack=True, dtype="str",skiprows=1)
  jtheta1, jtheta2, jbond1, jbond2  = np.loadtxt(dbdir+"cfangle.prm", usecols=(1, 2, 3, 4), unpack=True, dtype="float",skiprows=1)

  types = []; types_r = []
  jparams = []; jparams_r = []

  for t1, t2, t3 in zip(type1, type2, type3):
    types.append(t1 + "_" + t2 + "_" + t3)
    types_r.append(t3 + "_" + t2 + "_" + t1)

  for jt1, jt2, jb1, jb2 in zip(jtheta1, jtheta2, jbond1, jbond2):
    jparams.append(" ".join(["%10.5f"%jt1, "%10.5f"%jt2, "%10.5f"%jb1, "%10.5f"%jb2]))
    jparams_r.append(" ".join(["%10.5f"%jt2, "%10.5f"%jt1, "%10.5f"%jb2, "%10.5f"%jb1]))
  
  smartsCFangleDict = dict(zip(types, jparams))
  smartsCFangleDict_r = dict(zip(types, jparams_r))
  ttypes, stypes = np.loadtxt(f"{fname}.type", usecols=(1,2), unpack=True, dtype="str")
  tinker2smarts = dict(zip(ttypes, stypes))
  
  with open(tinkerkey, "a") as f:
    f.write("# cflux-a parameter assigned from database\n")
    for line in lines:
      if "angle " in line:
        dd = line.split()
        if (dd[1] in ttypes) and (dd[2] in ttypes) and (dd[3] in ttypes):
          s1 = tinker2smarts[dd[1]]
          s2 = tinker2smarts[dd[2]]
          s3 = tinker2smarts[dd[3]]
          comb1 = s1 + "_" + s2 + "_" + s3
          comb2 = s3 + "_" + s2 + "_" + s1
          if comb1 in smartsCFangleDict:
            f.write("cflux-a %s %s %s %s\n"%(dd[1], dd[2], dd[3], smartsCFangleDict[comb1]))
            print(GREEN + "CFlux parameters found for angle %s-%s-%s"%(dd[1], dd[2], dd[3]) + ENDC)
          elif comb2 in smartsCFangleDict:
            f.write("cflux-a %s %s %s %s\n"%(dd[1], dd[2], dd[3], smartsCFangleDict_r[comb2]))
            print(GREEN + "CFlux parameters found for angle %s-%s-%s"%(dd[1], dd[2], dd[3]) + ENDC)
          else:
            print(RED + "CFlux parameters NOT found for angle %s-%s-%s"%(dd[1], dd[2], dd[3]) + ENDC)
  return
 
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'xyz', help = "tinker xyz file", required=True)  
  parser.add_argument('-k', dest = 'key', help = "tinker prm file", required=True)  
  parser.add_argument('-p', dest = 'potent', help = "potential energy term", required=True, choices=["CF", "POLAR", "NONBONDED"])  
  args = vars(parser.parse_args())
  inp = args["xyz"]
  key = args["key"]
  potent = args["potent"]
  global dbdir
  dbdir="/home/liuchw/amoebaplus.data/"
  genAtomType(inp, potent)
  fname, _ = os.path.splitext(inp)
  if (potent.upper() == "POLAR"):
    assignPolar(fname, key)
  elif (potent.upper() == "CF"):
    assignCFlux(fname, key) 
  elif (potent.upper() == "NONBONDED"):
    assignNonbonded(fname, key) 
  else:
    print(RED + f"{potent} term not supported!" + ENDC)
  return

if len(sys.argv) == 1:
  print('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
else:
  main()
