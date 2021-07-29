
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import numpy as np

# color
RED = '\33[91m'
GREEN = '\33[92m'
ENDC = '\033[0m'

# read atoms and coordinate from txyz file
def readTXYZ(txyz):
  atoms  = np.loadtxt(txyz, usecols=(1,), dtype='str', unpack=True, skiprows=1)
  coords = np.loadtxt(txyz, usecols=(2,3,4), dtype='float', unpack=False, skiprows=1)
  return atoms,coords

# calculate distance between two atoms 
def distance(coord1, coord2):
  return np.sqrt(np.square(np.array(coord1)-np.array(coord2)).sum()) 

def writeCOM(atoms, coords, fname):
  with open(fname,"w") as f:
    nproc = "%Nproc=8\n" 
    memory = "%Mem=50GB\n"
    chk = "%Chk=" + "%s\n"%fname.replace(".com", ".chk")
    keywords = "#p MP2/aug-cc-pvtz SP Density=MP2 SCF=Save Charge NoSymm MaxDisk=100GB \n"
    comment = "SP job with an external charge\n"
    chgspin = f" {charge} 1  \n"
    f.write(chk + memory + nproc + keywords + "\n" + comment + "\n" + chgspin)
    for atom, coord in zip(atoms, coords):
      f.write("%3s %10.6f%10.6f%10.6f\n"%(atom, coord[0], coord[1], coord[2]))
  return

# get ESP values
def getQMESP(fname):
  QMfinished = False
  log = prefix + ".log"
  if os.path.isfile(log):
    for line in open(log).readlines():
      if "Normal termination" in line:
        QMfinished = True

  if (not os.path.isfile(f"{prefix}.fchk")) and QMfinished:
    cmdstr = f"{gaudir}/formchk {prefix}.chk"
    os.system(cmdstr)
  if (not os.path.isfile(f"{prefix}.cube")) and QMfinished:
    cmdstr = f"{gaudir}/cubegen 0 potential=MP2 {prefix}.fchk {prefix}.cube -5 h < {prefix}.grid"
    os.system(cmdstr)
  if (not os.path.isfile(f"{prefix}.pot")) and QMfinished:
    cmdstr = f"{tinkerdir}/potential.x 2 {prefix}.cube"
    os.system(cmdstr)
   
  # get pot for new grid from -esp.fchk

  name = fname + "-esp_" + prefix.split("_")[-1]
  cmdstr = f"{gaudir}/cubegen 0 potential=MP2 {fname}-esp.fchk {name}.cube -5 h < {prefix}.grid"
  os.system(cmdstr)
  
  cmdstr = f"{tinkerdir}/potential.x 2 {name}.cube "
  os.system(cmdstr)
  
  qmpol = np.loadtxt(prefix + ".pot", usecols=(-1), unpack=True, skiprows=1)
  qmref = np.loadtxt(name + ".pot", usecols=(-1), unpack=True, skiprows=1)
  np.savetxt(prefix + "_qm_polarization.pot", qmpol-qmref, fmt='%15.4f')
  return
    
def generate(d):
  atoms, coords = readTXYZ(txyz)
  vector = np.zeros(3)
  lines = open(txyz).readlines()

  if not perpend:
    for j in range(1,len(patom)):
      vector += (coords[patom[0]-1] - coords[patom[j]-1])
  else:
    if len(patom) == 3:
      v1 = coords[patom[0]-1] - coords[patom[1]-1]
      v2 = coords[patom[0]-1] - coords[patom[2]-1]
      vector = np.cross(v1,v2)
    else:
      sys.exit(RED + "Error: two atoms are needed to determine the direction of the probe" + ENDC)
  
  vector = vector/np.linalg.norm(vector)
  v = coords[patom[0]-1] + d*vector
  fname = prefix + ".com" 
  writeCOM(atoms, coords, fname)
  with open(fname, "a") as f:
    f.write("\n")
    f.write(" %10.6f%10.6f%10.6f%6.3f\n\n"%(v[0], v[1], v[2], 0.01))
    print(GREEN + "Generated %s"%fname + ENDC)
  fname = prefix + ".xyz"
  with open(fname, "w") as f:
    f.write("  %4s\n"%(len(atoms)+1))
    for line in lines[1:]:
      f.write(line)
    f.write("  %4s%3s  %12.6f%12.6f%12.6f      999\n"%(len(atoms)+1, "X", v[0], v[1], v[2]))
    print(GREEN + "Generated %s"%fname + ENDC)
  return

# get ESP values
def getMMESP(fname):
  if os.path.isfile("../final.key"):
    lines = open("../final.key").readlines()
    with open(f"{prefix}.key_0", 'w') as f:
      for line in lines:
        if 'polarize ' in line:
          d = line.split()
          line = '  '.join(d[0:2] + ['0.000'] + d[3:]) + '\n'
        else:
          line = line
        f.write(line)

      # write Probe
      f.write('atom 999 999 X "Probe Charge " 1 1.000 0\n')
      f.write('vdw 999 0.01 0.000\n')
      f.write('multipole 999 0.125\n')
      f.write('              0.000 0.000 0.000\n')
      f.write('              0.000\n')
      f.write('              0.000 0.000\n')
      f.write('              0.000 0.000 0.000\n')
  else:
    sys.exit("final.key does not exist")
  cmdstr = f"{tinkerdir}/potential.x 3 {prefix}.xyz -k {prefix}.key_0 Y > /dev/null "
  os.system(cmdstr)
  
  lines = open(f"{prefix}.pot").readlines()
  return 

# get ESP values
def getMMGRID(fname):
  if os.path.isfile("../final.key"):
    lines = open("../final.key").readlines()
    with open(f"{prefix}.key_0", 'w') as f:
      for line in lines:
        if 'polarize ' in line:
          d = line.split()
          line = '  '.join(d[0:2] + ['0.000'] + d[3:]) + '\n'
        else:
          line = line
        f.write(line)

      # write Probe
      f.write('atom 999 999 X "Probe Charge " 1 1.000 0\n')
      f.write('vdw 999 0.01 0.000\n')
      f.write('multipole 999 0.125\n')
      f.write('              0.000 0.000 0.000\n')
      f.write('              0.000\n')
      f.write('              0.000 0.000\n')
      f.write('              0.000 0.000 0.000\n')
  else:
    sys.exit("final.key does not exist")
  cmdstr = f"{tinkerdir}/potential.x 1 {prefix}.xyz -k {prefix}.key_0 "
  os.system(cmdstr)
  return 
  
def main():
  if len(sys.argv) == 1:
    sys.exit(RED + "run with -h option to see the usage"+ ENDC)
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'txyz', required=True, help="tinker xyz file not containing the probe ion")  
  parser.add_argument('-m', dest = 'mode', required=False, help="0: generating com/xyz files; 1:postprocessing", type=int, default=0)  
  parser.add_argument('-c', dest = 'charge', required=True, help="total charge of the system, not including probe charge")  
  parser.add_argument('-f', dest = 'prefix', required=True, help="prefix of generated filename")  
  parser.add_argument('-p', dest = "prob_atom", nargs='+', required=True, type=int, help="probe atom followed by directional atoms, e.g. 1 2 3: use bisector of 1-2 and 1-3") 
  parser.add_argument('-d', dest = 'distance', required=True, type=float, help="distance from the probing atom, e.g. 3 4 5: three distances at 3, 4 and 5 A")
  parser.add_argument('-a', dest = 'adaptingdistance', type=int, help="adapting distance for potential comparison purpose", default=0)
  args = vars(parser.parse_args())
  global gaudir,tinkerdir,txyz
  global prefix,charge,adapt
  global perpend,patom
  gaudir = "/opt/g09gh/gaussian/g09"
  tinkerdir = "$TINKER89"
  txyz = args["txyz"]
  mode  = args["mode"] 
  charge = args["charge"]
  prefix = args["prefix"]
  patom = args["prob_atom"]
  dist  = args["distance"]
  adapt  = args["adaptingdistance"]

  perpend = False
  for p in patom[1:]:
    if p < 0:
      perpend = True
      break
  
  fname = txyz.split(".")[0]
  # generating mode
  if mode == 0:
    generate(dist)
    getMMGRID(fname)
  
  # postprocessing mode
  elif mode == 1:
    getQMESP(fname)
  else:
    sys.exit(RED + "Error: only mode 0 and 1 are supported" + ENDC)
      
  return 

if __name__ == "__main__":
  main()
