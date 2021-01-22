
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import argparse
import os,sys,time
import subprocess
import numpy

# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i',     dest = 'input', required=True, help="input filename")  
  parser.add_argument('-o',     dest = 'output', default=None, help="output filename. Optional")
  parser.add_argument('-it',    dest = 'inType',  required=True, choices = ["xyz", "txyz", "g09", "qcout", "mol", "mol2", "sdf", "psiout"], help="input file type")
  parser.add_argument('-ot',    dest = 'outType', required=True, choices = ["xyz", "qcin", "psi4", "com", "txyz"], help="output file type")
  parser.add_argument('-q',     dest = 'QM', default="HF", type=str.upper, help="QM method")
  parser.add_argument('-b',     dest = 'basis',  default = "STO-3G", help="basis function for quantum job", type=str.lower)
  parser.add_argument('-c',     dest = 'charge', default = "0", help="total charge of the whole system", type=str.lower)
  parser.add_argument('-c1',    dest = 'charge1', default = None, help="total charge of the first molecule")
  parser.add_argument('-c2',    dest = 'charge2', default = None, help="total charge of the second molecule")
  parser.add_argument('-s',     dest = 'spin', default = "1", help="total spin of the whole system")
  parser.add_argument('-s1',    dest = 'spin1', default = None, help="total spin of the first molecule")
  parser.add_argument('-s2',    dest = 'spin2', default = None, help="total spin of the second molecule")
  parser.add_argument('-n1',    dest = 'number1', default = None, help="number of atoms of the first molecule")
  parser.add_argument('-bsse',  dest = 'bsse', default = None, help = "use bsse correction, e.g. CP")
  parser.add_argument('-j',     dest = 'jobType', default="sp", type=str.lower, help="job type, can be opt, sp, freq, cbs, sapt, opt+freq, dipole, esp. Default: sp")
  parser.add_argument('-d',     dest = 'disk', default = "10GB", help="size of disk to be used. Default: 10GB")
  parser.add_argument('-m',     dest = 'memory', default = "10GB", help="size of memory to be used. Default: 10GB")
  parser.add_argument('-n',     dest = 'nproc', default = "8", help="number of cores for Gaussian job. Default: 8")
  parser.add_argument('-conv',  dest = 'converge', default = " ", help="string, converge cariteria for gaussian opt job. Default: Gaussian default")
  parser.add_argument('-at',    dest = 'atomtype', nargs='+', default = None, help="atom types for txyz file. Please provide in a row")
  parser.add_argument('-tp',    dest = 'template', default = None, help="template txyz file used to convert to new txyz")
  args = vars(parser.parse_args())

  fi = args["input"]
  fo = args["output"]
  ti = args["inType"].upper()
  to = args["outType"].upper()
  if (fo == None):
    fo = os.path.splitext(fi)[0] + "." + to.lower()
    print(fo)
  cg = args["charge"]
  qm = args["QM"].upper()
  if qm == "CCSD_T":
    qm = "CCSD(T)"
  bf = args["basis"]
  sp = int(args["spin"])
  jt = args["jobType"]
  dk = args["disk"]
  me = args["memory"]
  np = args["nproc"]
  n1 = args["number1"]
  bsse = args["bsse"]
  c1 = args["charge1"]
  c2 = args["charge2"]
  s1 = args["spin1"]
  s2 = args["spin2"]
  at = args["atomtype"]
  tp = args["template"]
  con = args["converge"]

  def ToXYZ(fi,fo):
    if ti == "XYZ":cmdstr = "babel -ixyz %s -oxyz %s"%(fi, fo)
    elif ti == "TXYZ":
      lines = open(fi).readlines()
      n = len(lines[0].split())
      if n == 1:
        with open("tmp.t","w") as ftxyz:
          ftxyz.write(lines[0][:-1] + " Comments \n")
          for line in lines[1:]:
            ftxyz.write(line)
        cmdstr = "babel -itxyz %s -oxyz %s"%("tmp.t", fo)
      else:
        cmdstr = "babel -itxyz %s -oxyz %s"%(fi, fo)
    elif ti == "PSIOUT":
      psiout2xyz(fi,fo)
      cmdstr = "echo 'write xyz file from psi4 log!'"
    elif (ti in ["G09", "QCOUT", "MOL", "MOL2", "SDF"]):
      cmdstr = "babel -i%s %s -oxyz %s"%(ti.lower(), fi, fo)
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)
    return

  def ToTXYZ(fi,fo,at):
    if ti == "XYZ":cmdstr = "babel -ixyz %s -otxyz %s"%(fi, "tmp.txyz")
    elif ti == "PSIOUT":
      psiout2xyz(fi, "tmp.x")
      cmdstr = "babel -ixyz %s -otxyz %s"%("tmp.x", "tmp.txyz")
    elif (ti in ["TXYZ", "G09", "QCOUT", "MOL", "MOL2", "SDF"]):
      cmdstr = "babel -i%s %s -otxyz %s"%(ti.lower(), fi, "tmp.txyz")
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)
    atomtypes = []
    # read atom types from txyz files
    if at != None:
      for a in at:
        if not (a.endswith(".xyz") or a.endswith(".txyz")):
          types = numpy.loadtxt(a, usecols=(0,), dtype="str", unpack=True)
          atomtypes += list(types)
        else: # tinker-xyz
          types = numpy.loadtxt(a, usecols=(5,), dtype="str", unpack=True, skiprows=1)
          atomtypes += list(types)
    # read atom types and connections from txyz file
    atomtypes_connections = []
    if tp != None:
      lines = open(tp).readlines()
      for line in lines[1:]:
        dd = line.split()
        atomtypes_connections.append(dd[5:])
      
    lines = open("tmp.txyz").readlines()
    with open(fo, "w")  as f:
      f.write("%6s Generated by lconvert.py\n"%lines[0].split()[0])
      for i in range(1,len(lines)):
        dd = lines[i].split()
        newline = "%6s%3s%12.6f%12.6f%12.6f"%(dd[0], dd[1], float(dd[2]), float(dd[3]), float(dd[4]))
        if atomtypes != []:
          dd[5] = atomtypes[i-1]
          newline += "   ".join(["%4s"%i for i in dd[5:]])
        elif atomtypes_connections != []:
          newline = newline + " " +  "   ".join(atomtypes_connections[i-1])
        else:
          print("MM2 atom types are used!")
        f.write(newline + "\n")
    os.remove("tmp.txyz")
    return

  def ToCOM(fi,fo):
    if ti == "XYZ":
      cmdstr = "babel -ixyz %s -ogau %s"%(fi, fo)
    elif ti == "PSIOUT":
      psiout2xyz(fi, 'tmp.x')
      cmdstr = "babel -ixyz %s -ogau %s"%("tmp.x", fo)
    elif (ti in ["TXYZ", "COM", "G09", "QCOUT", "MOL", "MOL2", "SDF"]):
      cmdstr = "babel -i%s %s -ogau %s"%(ti.lower(), fi, fo)
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)

    fname, _ = os.path.splitext(fo)
    chk   = "%chk=" + fname + ".chk\n"
    nproc = "%Nproc="+np + "\n" 
    mem   = "%Mem="+me +"\n"
    
    if (jt.upper() == "OPT"):
      if con == ' ':
        extra =  "opt(calcFC, maxcycle=400) IOP(5/13=1) \n"
      else:
        extra =  "opt(calcFC, maxcycle=400, %s) IOP(5/13=1) \n"%con
    elif (jt.upper() == "OPT+FREQ"):
      extra = "IOP(5/13=1)\n"
    else:
      extra = '\n'

    if not bsse:
      counterpoise = ' '
    else:
      counterpoise = "counterpoise=2 " 

    if jt.upper() == "OPT+FREQ":
      key = " ".join(["#P", qm+"/"+bf, "opt", "freq", "MaxDisk=%s"%dk, extra])
    elif jt.upper() == "ESP":
      key = " ".join(["#P", qm+"/"+bf, " SP Density=Current SCF=Save NoSymm", "MaxDisk=%s"%dk, extra])
    elif jt.upper() == "OPT":
      key = " ".join(["#P", qm+"/"+bf,"NoSymm MaxDisk=%s"%dk, extra])
    else:
      key = " ".join(["#P", qm+"/"+bf, jt, "NoSymm MaxDisk=%s"%dk, extra])
    comment = " Generated by lconvert.py for %s job\n"%jt
    chgspin = " ".join([str(cg), str(sp), "\n"])
    
    with open(fo + "_2", 'w') as fout:
      fout.write(chk + mem + nproc + key + "\n" + comment + "\n" + chgspin)
      if not bsse:
        for line in open(fo).readlines()[5:]: 
          fout.write(line)
      else:
        lines = open(fo).readlines()[5:]
        for i in range(len(lines)-1):
          if i < int(n1):
            fout.write(lines[i].split("\n")[0] + "  1 \n")
          else:
            fout.write(lines[i].split("\n")[0] + "  2 \n")
        fout.write("\n")
    os.rename(fo+"_2", fo)
    return

  def ToQCHEM(fi,fo):
    if ti == "XYZ":
      cmdstr = "babel -ixyz %s -oqcin %s"%(fi, fo)
    elif ti == "PSIOUT":
      psiout2xyz(fi, "tmp.x")
      cmdstr = "babel -ixyz %s -oqcin %s"%("tmp.x", fo)
    elif (ti in ["TXYZ", "COM", "G09", "QCOUT", "MOL", "MOL2", "SDF"]):
      cmdstr = "babel -i%s %s -oqcin %s"%(ti.lower(), fi, fo)
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)

    basis   = " BASIS=%s\n"%bf
    jobtype = " JOB_TYPE=%s\n"%jt
    method  = " METHOD=%s\n"%qm
    chgspin = cg + " " + sp + "\n"
    with open(fo + "_2", 'w') as fout:
      fout.write("$comment\n Generated by lconvert.py for %s job\n$end\n\n$molecule\n"%jt)
      fout.write(chgspin)
      for line in open(fo).readlines()[6:-3]:
        fout.write(line)
      fout.write("$rem\n" + method + basis + jobtype + "$end\n")
    os.rename(fo+"_2", fo)
    return

  def ToPSI4(fi,fo):
    def XYZ2PSI4(xyz,psi4):
      atoms  = numpy.loadtxt(xyz, usecols=(0,), dtype='str', unpack=True, skiprows=2)
      coords = numpy.loadtxt(xyz, usecols=(1,2,3), dtype='float', unpack=True, skiprows=2)
      with open(psi4,'w') as fout:
        fout.write("#Generated by lconvert.py for %s job\n\n"%jt)
        if (bsse == None) and (jt.upper() != "SAPT"):
          chgspin = str(cg) + " " + str(sp)
          fout.write("memory %s %s\n\nmolecule  {\n%s\n"%(me[:-2], me[-2:], chgspin))
          for n in range(len(atoms)):
            fout.write("%3s   %12.6f%12.6f%12.6f\n"%(atoms[n],float(coords[0][n]),float(coords[1][n]),float(coords[2][n])))
        else: 
          chgspin1 = c1 + " " + s1
          chgspin2 = c2 + " " + s2
          fout.write("memory %s %s\n\nmolecule  {\n%s\n"%(me[:-2], me[-2:], chgspin1))
          for n in range(int(n1)):
            fout.write("%3s   %12.6f%12.6f%12.6f\n"%(atoms[n],float(coords[0][n]),float(coords[1][n]),float(coords[2][n])))
          fout.write("--\n%s\n"%(chgspin2))
          for n in range(int(n1), len(atoms)):
            fout.write("%3s   %12.6f%12.6f%12.6f\n"%(atoms[n],float(coords[0][n]),float(coords[1][n]),float(coords[2][n])))
        fout.write("units angstrom\nno_reorient\n")
        fout.write("symmetry c1\n}\n\n")
        fout.write("set {\nscf_type DF\n")
        if qm == "MP2":
          fout.write("mp2_type DF\ne_convergence 7\nreference rhf\n}\n\n")
          if jt.upper() == "CBS": 
            if not bsse:
              if bf.upper() == "CC-PVTZ":
                fout.write("energy('mp2/cc-pv[tq]z')\n")
              elif bf.upper() == "AUG-CC-PVTZ":
                fout.write("energy('mp2/aug-cc-pv[tq]z')\n")
              else:
                sys.exit("MP2 CBS with %s is not supported!"%bf.upper())
            else:
              if bf.upper() == "CC-PVTZ":
                fout.write("energy('mp2/cc-pv[tq]z', bsse_type='%s')\n"%bsse)
              elif bf.upper() == "AUG-CC-PVTZ":
                fout.write("energy('mp2/aug-cc-pv[tq]z', bsse_type='%s')\n"%bsse)
              else:
                sys.exit("MP2 CBS BSSE with %s is not supported!"%bf.upper())
          elif jt.upper() == "SP":
            if not bsse:
              fout.write("energy('mp2/%s')\n"%bf)
            else:
              fout.write("energy('mp2/%s', bsse_type='%s')\n"%(bf,bsse))
          elif jt.upper() == "OPT":
            fout.write("set opt_coordinates cartesian\n")
            fout.write("set g_convergence GAU\n")
            fout.write("set GEOM_MAXITER 400\n")
            fout.write("optimize('mp2/%s')\n"%bf)
            fout.write("\n")
          else:
            sys.exit("MP2 with %s is not supported!"%jt.upper())
        elif qm == "MP2D":
          fout.write("mp2_type DF\ne_convergence 7\nreference rhf\n}\n\n")
          if jt.upper() == "CBS": 
            if not bsse:
              if bf.upper() == "CC-PVTZ":
                fout.write("energy('mp2d/cc-pv[tq]z')\n")
              elif bf.upper() == "AUG-CC-PVTZ":
                fout.write("energy('mp2d/aug-cc-pv[tq]z')\n")
              else:
                sys.exit("MP2D CBS with %s is not supported!"%bf.upper())
            else:
              if bf.upper() == "CC-PVTZ":
                fout.write("energy('mp2d/cc-pv[tq]z', bsse_type='%s')\n"%bsse)
              elif bf.upper() == "AUG-CC-PVTZ":
                fout.write("energy('mp2d/aug-cc-pv[tq]z', bsse_type='%s')\n"%bsse)
              else:
                sys.exit("MP2D CBS BSSE with %s is not supported!"%bf.upper())
          elif jt.upper() == "SP":
            if not bsse:
              fout.write("energy('mp2d/%s')\n"%bf)
            else:
              fout.write("energy('mp2d/%s', bsse_type='%s')\n"%(bf,bsse))
          elif jt.upper() == "OPT":
            fout.write("optimize('mp2d/%s')\n"%bf)
          else:
            sys.exit("MP2D with %s is not supported!"%jt.upper())
        elif qm == "CCSD(T)":
          fout.write("mp2_type DF\ne_convergence 7\nreference rhf\n}\n\n")
          if jt.upper() == "CBS":
            if not bsse:
              fout.write("energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')\n")
            else:
              fout.write("energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', bsse_type='%s')\n"%bsse)
          if jt.upper() == "SP":
            if not bsse:
              fout.write("energy('ccsd(t)/%s')\n"%bf)
            else:
              fout.write("energy('ccsd(t)/%s', bsse_type='%s')\n"%(bf, bsse))
        else: 
          fout.write("e_convergence 7\nreference rhf\n\n\n")
          if jt.upper() == "SP":
            if not bsse:
              fout.write("energy('%s/%s')\n"%(qm,bf))
            else:
              fout.write("energy('%s/%s', bsse_type='%s')\n"%(qm,bf, bsse))

        if jt.upper() == "SAPT":
          fout.write("basis %s\n"%bf)
          fout.write("freeze_core True\n")
          fout.write("guess SAD\n}\n")
          fout.write("energy('sapt2+')\n")
        if jt.upper() == "DIPOLE":
          fout.write('set PROPERTIES_ORIGIN ["COM"]\n')
          fout.write("properties('PBE0/%s', properties=['dipole'], title='Acetate')\n"%bf)
      return

    rmstr = "rm -rf tmp.xyz"
    subprocess.run(rmstr, shell=True)
    if ti == "XYZ":
      cmdstr = "cp %s tmp.xyz"%fi
    elif ti == "PSIOUT":
      psiout2xyz(fi, "tmp.xyz")
      cmdstr = "echo 'write xyz file from psi4 log!'"
    elif (ti in ["TXYZ", "COM", "G09", "QCOUT", "MOL", "MOL2", "SDF"]):
      cmdstr = "babel -i%s %s -oxyz %s"%(ti.lower(), fi, "tmp.xyz")
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)

    if os.path.isfile("tmp.xyz"):
      XYZ2PSI4("tmp.xyz", fo)
    return
    
  def psiout2xyz(fin, fout):
    lines = open(fin).readlines()
    for n in range(len(lines)):
      if "Final optimized geometry " in lines[n]:break
    atoms = []
    coords = []
    for n in range(n+6, len(lines)):
      if lines[n] == "\n":break
      dd = lines[n].split()
      atoms.append(dd[0])
      coords.append(dd[1:4])
    with open(fout, "w") as fxyz:
      fxyz.write("%3s\nconverted from %s\n"%(len(atoms), fin))
      for n in range(len(atoms)):
        fxyz.write("%3s%12.6f%12.6f%12.6f\n"%(atoms[n], float(coords[n][0]), float(coords[n][1]), float(coords[n][2])))
    return
          
  if (to == "COM"):
    ToCOM(fi, fo)
  elif (to == "PSI4"):
    ToPSI4(fi, fo)
  elif (to == "QCIN"):
    ToQCHEM(fi, fo)
  elif (to == "XYZ"):
    ToXYZ(fi, fo)
  elif (to == "TXYZ"):
    ToTXYZ(fi, fo, at)
  else:
    print(RED + "File format '%s' not supported!"%to + ENDC)
  return

if __name__ == "__main__":
  if len(sys.argv) == 1:
    print(GREEN + " An example: " + ENDC + "python lconvert.py -i water.xyz -it xyz [-o water.com] -ot com -j opt -q MP2 -b cc-pvtz -c 0 -s 1 -d 10GB -m 20GB -n 10")
    sys.exit(RED + " For full usage, please run: " +ENDC + GREEN + "python lconvert.py -h" + ENDC)
  main()
