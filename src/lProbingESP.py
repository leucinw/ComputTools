
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import numpy as np

RED = '\33[91m'
GREEN = '\33[92m'
ENDC = '\033[0m'

def readTXYZ(txyz):
  atoms  = np.loadtxt(txyz, usecols=(1,), dtype='str', unpack=True, skiprows=1)
  coords = np.loadtxt(txyz, usecols=(2,3,4), dtype='float', unpack=False, skiprows=1)
  return atoms,coords

def distance(coord1, coord2):
  return np.sqrt(np.square(np.array(coord1)-np.array(coord2)).sum()) 

def psi4ESP():
  fname = f"{prefix}.psi4"
  atoms,coords = readTXYZ(txyz)
  with open(fname,"w") as f:
    f.write("import shutil\n\n")
    chgspin = f" {charge} 1  \n"
    f.write("molecule {\n")
    f.write(chgspin)
    for atom, coord in zip(atoms, coords):
      f.write("%3s %10.6f%10.6f%10.6f\n"%(atom, coord[0], coord[1], coord[2]))
    f.write("units angstrom\n")
    f.write("no_reorient\n")
    f.write("symmetry c1\n}\n\n")
    f.write("memory 48GB\n") 
    f.write("set_num_threads(8)\n")
    f.write("set maxiter 300\n")
    f.write("set freeze_core True\n")
    f.write('set PROPERTIES_ORIGIN ["COM"]\n\n')
    f.write("E, wfn = properties('MP2/aug-cc-pVTZ', properties=['GRID_ESP'], return_wfn=True)\n")
    f.write(f'fchk(wfn, "{prefix}.fchk")\n')
    f.write(f"wfn.to_file('{prefix}.npy')\n")
    f.write(f'shutil.move("grid_esp.dat", "{prefix}.grid_esp.dat")\n')
    f.write("clean()\n")
  return

def psi4ESP_prob(atoms, coords, ff, probcoord):
  fname = ff + ".psi4"
  with open(fname,"w") as f:
    f.write("import shutil\n\n")
    chgspin = f" {charge} 1"
    f.write("molecule {\n")
    f.write(chgspin + "\n")
    for atom, coord in zip(atoms, coords):
      f.write("%3s %10.6f%10.6f%10.6f\n"%(atom, coord[0], coord[1], coord[2]))
    f.write("units angstrom\n")
    f.write("no_reorient\n")
    f.write("symmetry c1\n}\n\n")
    f.write("memory 48GB\n") 
    f.write("set_num_threads(8)\n")
    f.write("set maxiter 300\n")
    f.write("set freeze_core True\n")
    f.write('set PROPERTIES_ORIGIN ["COM"]\n\n')
    f.write("Chrgfield = QMMM()\n")
    c = probcoord
    f.write("Chrgfield.extern.addCharge(0.125000,%10.6f,%10.6f,%10.6f)\n"%(c[0], c[1], c[2]))
    f.write("psi4.set_global_option_python('EXTERN', Chrgfield.extern)\n")
    f.write("E, wfn = properties('MP2/aug-cc-pVTZ', properties=['GRID_ESP'], return_wfn=True)\n")
    f.write(f'fchk(wfn, "{ff}.fchk")\n')
    f.write(f"wfn.to_file('{ff}.npy')\n\n")
    f.write(f'shutil.move("grid_esp.dat", "{ff}.grid_esp.dat")\n')
    f.write("clean()\n")
  return

def oneProbe(ff, probatoms):
  atoms, coords = readTXYZ(txyz)
  vector = np.zeros(3)
  lines = open(txyz).readlines()
  
  patom = []
  for p in probatoms.split():
    patom.append(int(p))

  perpend = False 
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
  v = coords[patom[0]-1] + dist*vector
  psi4ESP_prob(atoms, coords, ff, v)
  fname = ff + ".xyz"
  with open(fname, "w") as f:
    f.write("  %4s\n"%(len(atoms)+1))
    for line in lines[1:]:
      f.write(line)
    f.write("  %4s%3s  %12.6f%12.6f%12.6f      999\n"%(len(atoms)+1, "X", v[0], v[1], v[2]))
    print(GREEN + "Generated %s"%fname + ENDC)
  
  lines = open(f"{prefix}.key").readlines()
  fname = ff + ".key"
  with open(fname, "w") as f:
    for line in lines:
      f.write(line)
    f.write("# ESP probe parameters\n")
    f.write('atom          999    999   PC     "Probe Charge        "         1     1.000    0\n')
    f.write("vdw       999  0.0100   0.0000\n")
    f.write("multipole   999                         0.12500\n")
    f.write("                                        0.00000    0.00000    0.00000\n")
    f.write("                                        0.00000\n")
    f.write("                                        0.00000    0.00000\n")
    f.write("                                        0.00000    0.00000    0.00000\n")
    f.write("polarize           999          0.0000     0.0100\n")
    print(GREEN + "Generated %s"%fname + ENDC)
  return


def multipleProbes():
  root = os.getcwd()
  atoms, elements, atypes = np.loadtxt(f"{prefix}.xyz", usecols=(0,1,5), dtype='str', skiprows=1, unpack=True)
  atom_element_dict = dict(zip(atoms,elements))
  connections = []
  for line in open(f"{prefix}.xyz").readlines()[1:]:
    d = line.split()
    connections.append(d[6:])
  
  tmp = []
  probes = []
  for a,e,t,c in zip(atoms,elements,atypes,connections):
    patom = c[0] 
    for cc in c:
      if t not in tmp:
        tmp.append(t)
        if atom_element_dict[cc] != "H":
          patom = cc
        fname = f"{prefix}_{e}{a}"
        probatoms = f"{a} {patom}"
        oneProbe(fname,probatoms)
        os.system(f"mkdir {e}{a}")
        os.system(f"mv {fname}* ./{e}{a}")
        os.system(f"ln -s {os.getcwd()}/{prefix}.grid ./{e}{a}/grid.dat")
        probes.append(fname)
  
  with open("probes", "w") as f:
    for p in probes:
      f.write(p + "\n")
  return

def getGridESP():
  cmdstr = f"{tinkerdir}/potential.x 1 {prefix}.xyz -k {prefix}.key 2>/dev/null"
  os.system(cmdstr)
  cmdstr = f"{tinkerdir}/potential.x 3 {prefix}.xyz -k {prefix}.key Y  2>/dev/null"
  os.system(cmdstr)
  return 

def getPolarESP():
  lines = open("probes").readlines()
  currdir = os.getcwd()
  for line in lines: 
    pre = line.split("_")[0]
    d = line.split("_")[-1].split("\n")[0]
    ff = line.split("\n")[0]
    pr = np.loadtxt("%s/%s.grid_esp.dat"%(d,ff))
    rf = np.loadtxt(f"{pre}.grid_esp.dat")
    xs,ys,zs = np.loadtxt("grid.dat", usecols=(0,1,2), unpack=True)
    hartree2kcal = 627.509
    with open("qm_polarization.pot", "w") as f:
      f.write(f"{len(pr):>10d} qm_polarization\n")
      for i in range(len(pr)):
        f.write("%10s%14.6f%14.6f%14.6f%14.6f\n"%(i+1, xs[i], ys[i], zs[i], hartree2kcal*(pr[i]-rf[i])))
    os.system(f"mv qm_polarization.pot ./{d}")
  return

def main():
  if len(sys.argv) == 1:
    sys.exit(RED + "run with -h option to see the usage"+ ENDC)
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'txyz', required=True, help="tinker xyz file not containing the probe ion")  
  parser.add_argument('-k', dest = 'key', required=True, help="tinker keyfile without probe ion")  
  parser.add_argument('-m', dest = 'mode', required=False, help="0: generating com/xyz files; 1:postprocessing", type=int, default=0)  
  parser.add_argument('-c', dest = 'charge', required=True, help="total charge of the system, not including probe charge")  
  parser.add_argument('-f', dest = 'prefix', required=True, help="prefix of probe filename")  
  parser.add_argument('-d', dest = 'distance', required=True, type=float, help="distance from the probing atom")
  args = vars(parser.parse_args())
  global gaudir,tinkerdir,txyz,key
  global prefix,charge,dist
  tinkerdir = "$TINKER89"
  txyz = args["txyz"]
  key = args["key"]
  mode  = args["mode"] 
  charge = args["charge"]
  prefix = args["prefix"]
  dist  = args["distance"]

  # generating mode
  if mode == 0:
    fname = txyz.split(".")[0]
    if fname != prefix:
      os.system(f"ln -s {txyz} {prefix}.xyz")
    fname = key.split(".")[0]
    if fname != prefix:
      os.system(f"ln -s {key} {prefix}.key")
    if not os.path.isfile(f"{prefix}.psi4"):
      psi4ESP()
    if not os.path.isfile(f"{prefix}.pot"):
      getGridESP()
      os.system(f"ln -s {prefix}.grid grid.dat")
    multipleProbes()
  
  # postprocessing mode
  elif mode == 1:
    getPolarESP()
  else:
    sys.exit(RED + "Error: either 0 or 1" + ENDC)
      
  return 

if __name__ == "__main__":
  main()
