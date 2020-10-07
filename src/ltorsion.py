
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# 1. Rotate torsional bond and generate QM input
# 2. Fitting to QM obtaining torsional parameters 

# Usage: python ltorsion.py xyz key qmMethod torInterval fitnonezero mode 

helpinfo =\
''' 
   1.0 run ltorsion.py with Mode 1: generate torsional files, including QM and MM\n
   1.5 run QM jobs (for acceleration, e.g., use an external job submitting script, lsub.py)\n
   2.0 run ltorsion.py with Mode 2: fitting torsional parameters to the QM energy generated in 1.5\n
   3.0 run ltorsion.py with Mode 3: generate energy profile comparison\n
'''
 
import sys,os,subprocess,time,argparse
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-xyz', dest = 'inputxyz', required=True, type=str)  
parser.add_argument('-key', dest = 'inputkey', required=True, type=str)  
parser.add_argument('-optmethod',  dest = 'optmethod', default = "MP2", type=str)
parser.add_argument('-spmethod',  dest = 'spmethod', default = "MP2", type=str)
parser.add_argument('-optbasis',  dest = 'optbasis', default = "6-311++G(d,p)", type=str)
parser.add_argument('-spbasis',  dest = 'spbasis', default = "aug-cc-pvtz", type=str)
parser.add_argument('-charge', dest = 'charge', default="0", type=str)
parser.add_argument('-spin', dest = 'spin', default="1", type=str)
parser.add_argument('-mode', dest = 'mode', required=True, type=int)
parser.add_argument('-nstruct', dest = 'nstruct', default=12, type=int)
parser.add_argument('-interval', dest = 'interval', default=30, type=int)
parser.add_argument('-fitnonezero',  dest = 'fitnonezero', type=bool, default=True)
parser.add_argument('-restrain',  dest = 'restrain', default = 1)
args = vars(parser.parse_args())

inputxyz = args["inputxyz"] 
inputkey = args["inputkey"] 
optmethod = args["optmethod"] 
spmethod = args["spmethod"] 
optbasis = args["optbasis"] 
spbasis = args["spbasis"] 
charge = args["charge"] 
spin = args["spin"]
mode = args["mode"] 
nstruct = args["nstruct"] 
interval = args["interval"] 
fitnonezero = bool(args["fitnonezero"])
restrain = args["restrain"]
hartree2kcal = 627.5094740631 
degree2radian = np.pi/180.0
RED = '\33[91m'
GREEN = '\33[92m'
ENDC = '\033[0m'

class Torsion():
  # init
  def __init__(self, torstring, angle):
    if angle >= 0:
      self.fname = "TOR-%s-ANG_+%03d"%(torstring,angle)  
    else:
      self.fname = "TOR-%s-ANG_-%03d"%(torstring,abs(angle))
    self.com = self.fname + ".com"
    self.xyz = self.fname + ".xyz"
    self.key = self.fname + ".key"
    self.log = self.fname + ".log"
    self.out = self.fname + ".out"
    self.angle = angle
    self.torstr = " ".join(torstring.split("-"))


  # write gaussian input, .com file
  def writeQM(self):
    Atoms = np.loadtxt(inputxyz, usecols=(1,) , dtype="str", unpack=True, skiprows=1)
    Xs, Ys,Zs = np.loadtxt(inputxyz, usecols=(2,3,4) , dtype="float", unpack=True, skiprows=1)
    if not os.path.isfile(self.com):
      with open(self.com, 'w') as f:
        f.write("%chk="+"%s.chk\n"%self.fname)
        f.write("%Nproc=8\n")
        f.write("%Mem=100GB\n")
        f.write("#p %s/%s MaxDisk=100GB IOP(5/13=1) Opt=(ModRedundant, MaxCyc=200)\n\n"%(optmethod, optbasis))
        f.write("Restrained Optimization for torsion %s\n\n"%self.fname)
        f.write("%s %s\n"%(charge, spin))
        for atom, x, y, z in zip(Atoms, Xs, Ys, Zs):
          f.write("%3s%12.6f%12.6f%12.6f\n"%(atom, x,y,z))
        f.write("\n")
        f.write("%s =%6.1f B\n"%(self.torstr, self.angle))
        f.write("%s F\n\n"%(self.torstr))
        f.write("--Link1--\n")
        f.write("%chk="+"%s.chk\n"%self.fname)
        f.write("%Nproc=8\n")
        f.write("%Mem=100GB\n")
        f.write("#p %s/%s Geom=Checkpoint MaxDisk=100GB SP\n\n"%(spmethod,spbasis))
        f.write("Single-point energy for torsion %s\n\n"%self.fname)
        f.write("%s %s\n\n"%(charge, spin))

  # getQM energy from Normal terminated .log file 
  def getQM(self):
    if not os.path.isfile(self.log):
      print(RED + self.log + "does not exist!" + ENDC)
    else:
      lines = open(self.log).readlines()
      if not ("Normal termination" in lines[-1]):
        e = 0.0
      else:
        combinedString = ''.join([line[1:-1] for line in lines[-200:]])
        if spmethod.upper() == "MP2": 
          e = -float(combinedString.split("MP2=-")[1].split("\\")[0])*hartree2kcal
        else:
          e = -float(combinedString.split("HF=-")[1].split("\\")[0])*hartree2kcal
    return e

  # write MM .key
  def writeMM(self):
    atmnos, atmtyps = np.loadtxt(inputxyz, usecols=(0,5) , dtype="int", unpack=True, skiprows=1)
    AtomTypeDict = {}
    for atmno, atmtyp in zip(atmnos, atmtyps):
      if atmno not in AtomTypeDict:
        AtomTypeDict[atmno] = atmtyp
    subprocess.run("cp %s %s"%(inputxyz, self.xyz), shell=True)
    lines = open(inputkey).readlines()
    nprm = 0
    with open(self.key, "w") as f:
      for line in lines:
        if "torsion " not in line:
          f.write(line)
        else:
          ss = line.split()
          t1,t2,t3 = float(ss[5]), float(ss[8]), float(ss[11])
          if (fitnonezero == True):
            if t1 != 0:
              ss[5] = "PRM_%02d"%nprm
              nprm += 1
            if t2 != 0:
              ss[8] = "PRM_%02d"%nprm
              nprm += 1
            if t3 != 0:
              ss[11] = "PRM_%02d"%nprm
              nprm += 1
          else:
            ss[5] = "PRM_%02d"%nprm
            ss[8] = "PRM_%02d"%(nprm+1)
            ss[11] = "PRM_%02d"%(nprm+2)
            nprm += 3
          f.write("   ".join(ss) + "\n")
    with open(self.key, "a") as f:
      f.write("\nrestrain-torsion %s 1.0 %6.1f\n"%(self.torstr, self.angle))
    return 
  # getMM energy from tinker .out file 
  def getMM(self):
    for line in open(self.out).readlines():
      if "Final Function Value :" in line:
        e = float(line.split()[-1])
    return e

  # get torsion list
  @staticmethod
  def getTorlist():
    atmnos, atmtyps = np.loadtxt(inputxyz, usecols=(0,5) , dtype="int", unpack=True, skiprows=1)
    AtomTypeDict = {}
    for atmno, atmtyp in zip(atmnos, atmtyps):
      if atmno not in AtomTypeDict:
        AtomTypeDict[atmno] = atmtyp
    if not os.path.isfile("ttt.out"):
      subprocess.run("analyze.x %s -k %s EP > ttt.out && wait"%(inputxyz,inputkey), shell=True)
    idx = 0
    lines = open("ttt.out").readlines()
    for line in lines: 
      if "Torsional Angle Parameters :" in line:
        idx = lines.index(line) + 4
    uniqueTorlistType = []
    uniqueTorlistAtom = []
    torPrmListAtom = []
    for i in range(idx, len(lines), 1):
      ss = lines[i].split()
      if len(ss) == 0:break
      else:
        torprm = [0.0, 0.0, 0.0]
        if "0/1" in ss:
          torprmidx = ss.index("0/1") - 1
          torprm[0] = float(ss[torprmidx])
        if "180/2" in ss:
          torprmidx = ss.index("180/2") - 1
          torprm[1] = float(ss[torprmidx])
        if "0/3" in ss:
          torprmidx = ss.index("0/3") - 1
          torprm[2] = float(ss[torprmidx])
        torlist0 = [int(ss[1]), int(ss[2]), int(ss[3]), int(ss[4])]
        torlist1 = "-".join(list(map(str,[AtomTypeDict[int(ss[1])], AtomTypeDict[int(ss[2])], AtomTypeDict[int(ss[3])], AtomTypeDict[int(ss[4])]])))
        torlist2 = "-".join(list(map(str,[AtomTypeDict[int(ss[4])], AtomTypeDict[int(ss[3])], AtomTypeDict[int(ss[2])], AtomTypeDict[int(ss[1])]])))
        if (torlist1 not in uniqueTorlistType) and (torlist2 not in uniqueTorlistType):
          uniqueTorlistType.append(torlist1)
          uniqueTorlistAtom.append(torlist0)
          torPrmListAtom.append(torprm)
    return uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom
  
  # get the initial torsion value
  @staticmethod
  def dihedral(n0, n1, n2, n3):
    Xs, Ys,Zs = np.loadtxt(inputxyz, usecols=(2,3,4) , dtype="float", unpack=True, skiprows=1)
    p0 = np.array([Xs[n0-1], Ys[n0-1], Zs[n0-1]]) 
    p1 = np.array([Xs[n1-1], Ys[n1-1], Zs[n1-1]]) 
    p2 = np.array([Xs[n2-1], Ys[n2-1], Zs[n2-1]]) 
    p3 = np.array([Xs[n3-1], Ys[n3-1], Zs[n3-1]]) 
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x)) 

  # get the initial torsion parameters from ttt.key 
  @staticmethod
  def refparams():
    p0 = []
    lines = open(inputkey).readlines()
    for line in lines:
      if "torsion " in line:
        ss = line.split()
        t1,t2,t3 = float(ss[5]), float(ss[8]), float(ss[11])
        if (fitnonezero == True):
          if t1 != 0:
            p0.append(t1) 
          if t2 != 0:
            p0.append(t2) 
          if t3 != 0:
            p0.append(t3) 
        else:
          p0 = p0 + [t1, t2, t3]
    return np.array(p0) 

  # write runMin.sh
  @staticmethod
  def writeRunMin():
    files = os.listdir()
    i = 0
    with open("runMin.sh", "w") as fo:
      for f in files:
        if ("TOR" in f) and (".xyz" in f):
          if i == 12:
            fo.write("minimize.x %s -k %s_2 0.01 > %s \n"%(f, f.replace(".xyz",".key"), f.replace(".xyz", ".out")))
            i = 0
          else:
            fo.write("minimize.x %s -k %s_2 0.01 > %s &\n"%(f, f.replace(".xyz",".key"), f.replace(".xyz", ".out")))
            i += 1
    return

# get QM input and runMin.sh
def DataPrep():
  uniqueTorlistAtom, _, _, = Torsion.getTorlist()
  subprocess.run("rm -f runMin.sh", shell=True)
  for tor in uniqueTorlistAtom:
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3]))
    start = int(nstruct*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <= -180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      tor.writeQM()
      tor.writeMM()
  Torsion.writeRunMin()
  return 

# cost Function 
def costFuncTor(params):
  subprocess.run("rm -f TOR*.xyz_2 && wait",shell=True)
  uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom = Torsion.getTorlist()
  files = os.listdir(os.getcwd())
  for f in files:
    if ("TOR" in f) and f.endswith(".key"):
      with open(f+"_2", "w") as key2:
        key2.write("openmp-threads 1\n")
        for line in open(f).readlines():
          if ("torsion " not in line) or ("restrain-torsion" in line):
            key2.write(line)
          else:
            sss = line.split()
            if "PRM_" in sss[5]:
              idx = int(sss[5].split("_")[1])
              prmstr = "PRM_%02d"%idx
              line = line.replace(prmstr, str("%15.8f"%params[idx]))
            if "PRM_" in sss[8]:
              idx = int(sss[8].split("_")[1])
              prmstr = "PRM_%02d"%idx
              line = line.replace(prmstr, str("%15.8f"%params[idx]))
            if "PRM_" in sss[11]:
              idx = int(sss[11].split("_")[1])
              prmstr = "PRM_%02d"%idx
              line = line.replace(prmstr, str("%15.8f"%params[idx]))
            key2.write(line)
  subprocess.run("sh runMin.sh",shell=True)
  grepstr = "grep 'Final Function Value :' TOR*.out > grep.dat" 
  nline1 = len(open("runMin.sh").readlines())
  while True:
    subprocess.run(grepstr, shell=True)
    nline2 = len(open("grep.dat").readlines())
    if nline1 == nline2:break
    else:
      time.sleep(1.0)
  QM = []; MM = []
  for tor,torT,prm in zip(uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom):
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3]))
    start = int(nstruct*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <=-180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      qm = tor.getQM()
      if (qm != 0.0):
        QM.append(qm)
        MM.append(tor.getMM())
  QM = np.array(QM) - np.array(QM).min() 
  MM = np.array(MM) - np.array(MM).min() 
  for i in range(len(QM)):
    if (QM[i] > 20.0) or (MM[i]>20.0):
      QM[i] = 0.0 
      MM[i] = 0.0
  print(GREEN + " current RMSE is %10.5f"%np.sqrt(np.square(QM-MM).mean()) + ENDC)
  return QM-MM

# fitting all torsion parameters together
def fittingTorsion():
  x0 = Torsion.refparams()
  lowerBound = []
  upperBound = []
  for x in x0:
    if abs(x) > 0:
      lowerBound.append(x-abs(x)*0.2)
      upperBound.append(x+abs(x)*0.2)
    else:
      lowerBound.append(-0.1)
      upperBound.append(0.1)
  if restrain == 1:
    print(RED + "Doing torsion fitting with restraint" + ENDC)
    ret = least_squares(costFuncTor, x0, bounds = (lowerBound, upperBound), verbose=2, diff_step=0.0001, ftol=0.0001, gtol=0.0001, xtol=0.0001)
  else:
    print(RED + "Doing torsion fitting without restraint" + ENDC)
    ret = least_squares(costFuncTor, x0, verbose=2, diff_step=0.0001, ftol=0.0001, gtol=0.0001, xtol=0.0001)

# plot torsion energy profile and write tinker.key
def plotData():
  uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom = Torsion.getTorlist()
  for tor,torT,prm in zip(uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom):
    QM = []; MM = []; Ang = []
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3]))
    start = int(nstruct*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <= -180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      qm = tor.getQM()
      if (qm != 0.0):
        QM.append(qm)
        MM.append(tor.getMM())
        Ang.append(tor.angle)
    QM = np.array(QM) - np.array(QM).min() 
    MM = np.array(MM) - np.array(MM).min()
    for i in range(len(QM)):
      if (QM[i] > 20.0) or (MM[i] >20.0):
        QM[i] = 0.0
        MM[i] = 0.0
    fig, ax = plt.subplots()

    QM = [x for _,x in sorted(zip(Ang,QM))]
    MM = [x for _,x in sorted(zip(Ang,MM))]
    Ang = sorted(Ang)
 
    ax.plot(Ang, QM, label="QM")
    ax.plot(Ang, MM, label="MM")
    ax.set(xlabel='Torsional Angle (deg.)', ylabel='Energy (kcal/mol)', title='Torsion-%s'%torstring)
    ax.grid()
    ax.legend(loc="best")
    fig.savefig("torsion_fitting-%s.png"%torstring)
    plt.close() 

  files = os.listdir(os.getcwd())
  for f in files:
    if ("TOR" in f) and f.endswith(".key_2"):break
  with open("tinker.key", "w") as ftinker:
    for line in open(f).readlines():
      if "restrain-torsion" not in line:
        ftinker.write(line)
  return
# Execute different MODE
if mode == 1:
  DataPrep()
# After mode 1, QM jobs should be run
# python ~/bin/lsub.py -i TOR*.com
elif mode == 2:
  fittingTorsion()
elif mode == 3:
  plotData()
else:
  print(GREEN + helpinfo + ENDC)
