
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import argparse
import os,sys,time
import subprocess
from colorama import Fore

def main():
  #===>>>
  parser = argparse.ArgumentParser()
  parser.add_argument('-i',  dest = 'input', required=True)  
  parser.add_argument('-o',  dest = 'output', required=True)
  parser.add_argument('-it', dest = 'inType', required=True, choices=["xyz", "txyz", "g09", "qcout", "mol", "mol2", "sdf"])
  parser.add_argument('-ot', dest = 'outType', required=True, choices = ["xyz", "qcin", "psi4", "com"])
  parser.add_argument('-q',  dest = 'QM', choices = ["HF", "MP2", "B3LYP", "wB97XD", "CCSD(T)"], default="HF")
  parser.add_argument('-b',  dest = 'basis',  default = "STO-3G")
  parser.add_argument('-c',  dest = 'charge', default = "0")
  parser.add_argument('-s',  dest = 'spin', default = "1")
  parser.add_argument('-j',  dest = 'jobType', choices = ["opt", "sp", "freq", "cbs", "sapt", "opt+freq"], default="sp")
  parser.add_argument('-d',  dest = 'disk', default = "10GB")
  parser.add_argument('-m',  dest = 'memory', default = "10GB")
  parser.add_argument('-n',  dest = 'nproc', default = "8")
  args = vars(parser.parse_args())

  #!==>
  #==========================================================================
  def ToXYZ(fi,fo):
    if ti == "XYZ":
      cmdstr = "babel -ixyz %s -oxyz %s"%(fi, fo)
    elif ti == "TXYZ":
      cmdstr = "babel -itxyz %s -oxyz %s"%(fi, fo)
    elif ti == "G09":
      cmdstr = "babel -ig09 %s -oxyz %s"%(fi, fo)
    elif ti == "QCOUT":
      cmdstr = "babel -iqcout %s -oxyz %s"%(fi, fo)
    elif ti == "MOL":
      cmdstr = "babel -imol %s -oxyz %s"%(fi, fo)
    elif ti == "MOL2":
      cmdstr = "babel -imol2 %s -oxyz %s"%(fi, fo)
    elif ti == "SDF":
      cmdstr = "babel -isdf %s -oxyz %s"%(fi, fo)
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)
    return
  #==========================================================================
  #!<==

  #!==>
  #==========================================================================
  def ToCOM(fi,fo):
    if ti == "XYZ":
      cmdstr = "babel -ixyz %s -ogau %s"%(fi, fo)
    elif ti == "TXYZ":
      cmdstr = "babel -itxyz %s -ogau %s"%(fi, fo)
    elif ti == "COM":
      cmdstr = "babel -icom %s -ogau %s"%(fi, fo)
    elif ti == "G09":
      cmdstr = "babel -ig09 %s -ogau %s"%(fi, fo)
    elif ti == "QCOUT":
      cmdstr = "babel -iqcout %s -ogau %s"%(fi, fo)
    elif ti == "MOL":
      cmdstr = "babel -imol %s -ogau %s"%(fi, fo)
    elif ti == "MOL2":
      cmdstr = "babel -imol2 %s -ogau %s"%(fi, fo)
    elif ti == "SDF":
      cmdstr = "babel -isdf %s -ogau %s"%(fi, fo)
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)

    fname, _ = os.path.splitext(fo)
    chk   = "%chk=" + fname + ".chk\n"
    nproc = "%Nproc="+np + "\n" 
    mem   = "%Mem="+me +"\n"
    
    if (jt.upper() == "OPT"):
      extra = "IOP(5/13=1) CalcFC\n"
    elif (jt.upper() == "OPT+FREQ"):
      extra = "IOP(5/13=1)\n"
    else:
      extra = '\n'

    if jt.upper() == "OPT+FREQ":
      key = " ".join(["#P", qm+"/"+bf, "opt", "freq", "MaxDisk=%s"%dk, extra])
    else:
      key = " ".join(["#P", qm+"/"+bf, jt, "MaxDisk=%s"%dk, extra])
    comment = " Generated by lconvert.py for %s job\n"%jt
    chgspin = " ".join([cg, sp, "\n"])
    
    with open(fo + "_2", 'w') as fout:
      fout.write(chk + mem + nproc + key + "\n" + comment + "\n" + chgspin)
      for line in open(fo).readlines()[5:]: 
        fout.write(line)
    os.rename(fo+"_2", fo)
    return
  #==========================================================================
  #!<==


  #!==>
  #==========================================================================
  def ToQCHEM(fi,fo):
    if ti == "XYZ":
      cmdstr = "babel -ixyz %s -oqcin %s"%(fi, fo)
    elif ti == "TXYZ":
      cmdstr = "babel -itxyz %s -oqcin %s"%(fi, fo)
    elif ti == "COM":
      cmdstr = "babel -icom %s -oqcin %s"%(fi, fo)
    elif ti == "G09":
      cmdstr = "babel -ig09 %s -oqcin %s"%(fi, fo)
    elif ti == "QCOUT":
      cmdstr = "babel -iqcout %s -oqcin %s"%(fi, fo)
    elif ti == "MOL":
      cmdstr = "babel -imol %s -oqcin %s"%(fi, fo)
    elif ti == "MOL2":
      cmdstr = "babel -imol2 %s -oqcin %s"%(fi, fo)
    elif ti == "SDF":
      cmdstr = "babel -isdf %s -oqcin %s"%(fi, fo)
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
  #==========================================================================
  #!<==

  #!==>
  #==========================================================================
  def ToPSI4(fi,fo):
    import numpy as np
    def XYZ2PSI4(xyz,psi4):
      atoms  = np.loadtxt(xyz, usecols=(0,), dtype='str', unpack=True, skiprows=2)
      coords = np.loadtxt(xyz, usecols=(1,2,3), dtype='float', unpack=True, skiprows=2)
      with open(psi4,'w') as fout:
        chgspin = cg + " " + sp 
        fout.write("#Generated by lconvert.py for %s job\n\n"%jt)
        fout.write("memory %s %s\n\nmolecule  {\n%s\n"%(me[:-2], me[-2:], chgspin))
        for n in range(len(atoms)):
          fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coords[0][n]),float(coords[1][n]),float(coords[2][n])))
        fout.write("units angstrom\nno_reorient\n")
        fout.write("symmetry c1\n}\n\n")
        fout.write("set {\nscf_type DF\n")
        if qm.upper() == "MP2":
          fout.write("mp2_type DF\ne_convergence 7\nreference rhf\n}\n\n")
          if jt.upper() == "CBS" and bf.upper()=="AUG-CC-PVTZ":
            fout.write("energy('mp2/aug-cc-pv[tq]z')\n")
          if jt.upper() == "CBS" and bf.upper()=="CC-PVTZ":
            fout.write("energy('mp2/cc-pv[tq]z')\n")
          if jt.upper() == "SP":
            fout.write("energy('mp2/%s')\n"%bf)
        if qm.upper() == "CCSD(T)":
          if jt.upper() == "CBS":
            fout.write("energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')\n")
          if jt.upper() == "SP":
            fout.write("energy('ccsd(t)/%s')\n"%bf)
      return

    rmstr = "rm -rf tmp.xyz"
    subprocess.run(rmstr, shell=True)
    if ti == "XYZ":
      cmdstr = "cp %s tmp.xyz"%ti
    elif ti == "TXYZ":
      cmdstr = "babel -itxyz %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "COM":
      cmdstr = "babel -icom %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "G09":
      cmdstr = "babel -ig09 %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "QCOUT":
      cmdstr = "babel -iqcout %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "MOL":
      cmdstr = "babel -imol %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "MOL2":
      cmdstr = "babel -imol2 %s -oxyz %s"%(fi, "tmp.xyz")
    elif ti == "SDF":
      cmdstr = "babel -isdf %s -oxyz %s"%(fi, "tmp.xyz")
    else:
      cmdstr = ("echo 'File format %s not supported!'"%ti)
    subprocess.run(cmdstr,shell=True)

    if os.path.isfile("tmp.xyz"):
      XYZ2PSI4("tmp.xyz", fo)
    return
  #==========================================================================
  #!<==
  
  #===>>> 
  fi = args["input"]
  fo = args["output"]
  ti = args["inType"].upper()
  to = args["outType"].upper()
  cg = args["charge"]
  qm = args["QM"]
  bf = args["basis"]
  sp = args["spin"]
  jt = args["jobType"]
  dk = args["disk"]
  me = args["memory"]
  np = args["nproc"]

  if (to == "COM"):
    ToCOM(fi, fo)
  elif (to == "PSI4"):
    ToPSI4(fi, fo)
  elif (to == "QCIN"):
    ToQCHEM(fi, fo)
  elif (to == "XYZ"):
    ToXYZ(fi, fo)
  else:
    print(Fore.RED + "File format '%s' not supported!"%to)
  return

if __name__ == "__main__":
  main()

