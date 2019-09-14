#!/usr/bin/env python

#  TXYZ2COM:  From TINKER xyz file to Gaussian input file
#  LOG2COM:   From Gaussian output file to Gaussian input file
#  COM2TXYZ:  From Gaussian input to TINKER XYZ file
#  TXYZ2XYZ:  From TINKER xyz file to common XYZ file (most software readable) 
#  TXYZ2TXYZ: From TINKER xyz file to a txyz file that all atoms are in order 
#  TXYZ2XYZi: From TINKER xyz file to common XYZ file without metal ion
#  COM2XYZ:   From Gaussian input file to common XYZ file
#  XYZ2QCHEM: From common XYZ file to Qchem file
#  COM2COM:   From Gaussian input file to new input file
#  ARC2ARC:   From TINKER traj file to traj file 
#  XYZ2TXYZ:  From common XYZ to Tinker XYZ file 
#  XYZ2COM:   From common xyz file to Gaussian input file
#  LOG2XYZ:   From Gaussian output file to normal xyz file 
#  XYZ2XYZ:   From normal xyz file to normal xyz file in different distances
#  XYZ2PSI4:  From normal xyz file to psi4 file(SAPT calculation) 

import sys
import os,string
from readChemFile import *
from operation import *

#
# XXX --> COM 
#
# 1. com to com
def COM2COM(COM):
  data = readCOM(COM)
  atoms = data[0];coord = data[1]
  fileName = COM.split(".com")[0]
  fout=open(fileName+".gjf",'w')
  fout.write("%chk="+"%s.chk\n"%fileName)
  fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(len(atoms)):
    fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

# 2. common xyz to com
def XYZ2COM(XYZ, optCartesian=False):
  QM='HF' 
  basis ='3-21G ' 
  nproc     ='4' 
  chrgMult='0 1'
  jobType = "Opt IOP(5/13=1)"
  extraKeyword = " NoSymm Guess=INDO MaxDisk=100GB "
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  fout=open(fileName+".com",'w')
  fout.write("%chk="+"%s.chk\n"%fileName)
  fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(len(atoms)):
    fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

def XYZ2COM_BSSE(XYZ, nfirst, chrgMult, optCartesian=False):
  QM='MP2' 
  basis ='cc-pvtz ' 
  nproc     ='8' 
  jobType = "Opt=tight "
  if optCartesian == True:
    jobType = "Opt=(tight,Cartesian) "
  extraKeyword = " counterpoise=2 IOP(5/13=1) scf=tight Symmetry=None MaxDisk=100GB "
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  fout=open(fileName+".com",'w')
  fout.write("%NoSave\n")
  fout.write("%chk="+"%s.chk\n"%fileName)
  fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(nfirst):
    fout.write("%3s(Fragment=1)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  for n in range(nfirst,len(atoms)):
    fout.write("%3s(Fragment=2)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

def XYZ2PSI4(XYZ, nfirst, chrgMult_1, chrgMult_2, basis):
  data = readXYZ(XYZ)
  atoms = data[0]
  coord = data[1]
  fout=open(XYZ.split(".xyz")[0]+".psi4","w")
  fout.write("#psi4\n#\nmemory 32 gb\n\nmolecule  {\n%s\n"%chrgMult_1)
  for n in range(nfirst):
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("--\n%s\n"%chrgMult_2)
  for n in range(nfirst, len(coord)):
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("units angstrom\nno_reorient  #important for SAPT in psi4, default\n")
  fout.write("symmetry c1  #important for SAPT in psi4, default\n}\n\n")
  if basis == "TZ":
    fout.write("set {\nbasis aug-cc-pVTZ\nscf_type DF\nfreeze_core True\n}\nenergy('sapt2+')\n")
  if basis == "DZ":
    fout.write("set {\nbasis aug-cc-pVDZ\nscf_type DF\nfreeze_core True\n}\nenergy('sapt2+')\n")
  return True

# 3. tinker xyz to com
def TXYZ2COM(TXYZ, nfirst, counterpoise=False):
  # scf=tight opt=tight counterpoise=2
  QM='MP2' 
  basis ='cc-pvtz ' 
  nproc     ='4' 
  chrgMult='0 1'
  jobType = "Opt=tight "
  extraKeyword = " scf=tight counterpoise=2 "
  
  #extraKeyword =' scf(fermi, tight) density(current) ' 
  #jobType = "polar "
  #jobType ='Opt(modred,gdiis)'
  #extraKeyword =' IOP(5/13=1) Int(grid=fine) SCRF(PCM,Solvent=Chloroform,Read) Test  ' 
  data = readTXYZ(TXYZ)
  atoms = data[0];coord = data[1]
  fileName = TXYZ.split(".txyz")[0]
  fout=open(fileName+".com",'w')
  fout.write("%chk="+"%s.chk\n"%fileName)
  if counterpoise==False:
    fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
    for n in range(len(atoms)):
      fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  else:
    extraKeyword = 'counterpoise=2 '
    chrgMult = '0,1 0,1 0,1 '
    fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
    for n in range(nfirst):
      fout.write("%3s(Fragment=1)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
    for n in range(nfirst,len(atoms)):
      fout.write("%3s(Fragment=2)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

# 4. gaussian log to com
def LOG2COM(LOG, counterpoise=False,optCartesian=False):
  pTable = {'1':'H','2':'He','3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O',\
                 '9':'F','10':'Ne','11':'Na','12':'Mg','13':'Al','14':'Si','15':'P',\
                 '16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca','22':'Ti','26':'Fe',\
                 '27':'Co','28':'Ni','29':'Cu','30':'Zn','35':'Br'}
  data = readLOG(LOG)
  if optCartesian == True:
    jobType = "Opt=(tight,Cartesian) "
  coord = data[0] 
  charge = data[1]
  multip = data[2] 
  fileName = LOG.split(".out")[0]
  fout=open(fileName+".com_1",'w')
  if counterpoise==False:
    fout.write("%chk="+"%s.chk\n"%(fileName))
    fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%3s%3s\n"%(jobType,jobType,charge,multip))
    for n in range(len(coord)):
      line=coord[n].split()
      fout.write("%3s             %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
  else:
    extraKeyword = 'counterpoise=2 '
    chrgMult = '0,1 0,1 0,1 '
    nfirst = len(coord)/2 # homo dimer 
    fout.write("%chk="+"%s.chk\n"%(fileName))
    fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
    for n in range(nfirst):
      line=coord[n].split()
      fout.write("%3s(Fragment=1)             %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
    for n in range(nfirst,len(coord)):
      line=coord[n].split()
      fout.write("%3s(Fragment=2)             %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
  fout.write("\n")
  fout.close()
  return True

def LOG2COM_BSSE(LOG, nfirst, optCartesian=False):
  pTable = {'1':'H','2':'He','3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O',\
                 '9':'F','10':'Ne','11':'Na','12':'Mg','13':'Al','14':'Si','15':'P',\
                 '16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca','22':'Ti','26':'Fe',\
                 '27':'Co','28':'Ni','29':'Cu','30':'Zn','35':'Br'}
  data = readLOG(LOG)
  coord = data[0]
  QM='MP2' 
  basis ='cc-pvtz ' 
  nproc     ='8' 
  jobType = "Opt=tight "
  if optCartesian == True:
    jobType = "Opt=(tight,Cartesian) "
  extraKeyword = " counterpoise=2 IOP(5/13=1) scf=tight Symmetry=None MaxDisk=100GB "
 
  fileName = LOG.split(".out")[0]
  fout=open(fileName+".com_1",'w')
  chrgMult = '0,1 0,1 0,1 '
  fout.write("%chk="+"%s.chk\n"%(fileName))
  fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(nfirst):
    line=coord[n].split()
    fout.write("%3s(Fragment=1)             %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
  for n in range(nfirst,len(coord)):
    line=coord[n].split()
    fout.write("%3s(Fragment=2)             %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
  fout.write("\n")
  fout.close()
  return True

def LOG2XYZ(LOG):
  pTable = {'1':'H','2':'He','3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O',\
                 '9':'F','10':'Ne','11':'Na','12':'Mg','13':'Al','14':'Si','15':'P',\
                 '16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca','22':'Ti','26':'Fe',\
                 '27':'Co','28':'Ni','29':'Cu','30':'Zn','35':'Br'}
  data = readLOG(LOG)
  coord = data[0]
  fileName = LOG.split(".out")[0]
  fout=open(fileName+".xyz",'w')
  fout.write("%3s\n#Generated by chemFileConvert.py\n"%len(coord))
  for n in range(len(coord)):
    line=coord[n].split()
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(pTable[line[1]],float(line[3]),float(line[4]),float(line[5])))
  fout.close()
  return True


def COM2XYZ(COM):
  data = readCOM(COM)
  atoms = data[0]
  coord = data[1]
  oFile = open(COM.split(".")[0]+".xyz","w")
  oFile.write("%3s\n#Generated by chemFileConvert.py\n"%len(atoms))
  for i in range(len(atoms)):
    oFile.write("%3s%14.7f%14.7f%14.7f\n"%(atoms[i], coord[i][0], coord[i][1], coord[i][2]))
  oFile.close()
  return True

def XYZ2XYZ(XYZ, nfirst, refFirst, refSecond, dist_ratio):
  data = readXYZ(XYZ)
  atoms = data[0]
  coord = data[1]
  coord_1 = coord[:nfirst]
  coord_2 = coord[nfirst:]
  coord = f_translate(coord_1, coord_2, refFirst, refSecond, dist_ratio)
  ratio = str("%.02f"%dist_ratio)
  print("You are doing distance "+ratio)
  oFile = open(XYZ.split(".")[0]+"_%s.xyz"%ratio,"w")
  oFile.write("%3s\n#Generated by chemFileConvert.py\n"%len(atoms))
  for i in range(len(atoms)):
    oFile.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[i], coord[i][0], coord[i][1], coord[i][2]))
  oFile.close()
  return True

def TXYZ2XYZ(TXYZ):
  data = readTXYZ(TXYZ)
  atoms = data[0];coord = data[1]
  fileName = TXYZ.split(".txyz")[0]
  fout=open(fileName+".xyz",'w')
  fout.write("%3s\n#Generated by chemFileConvert.py\n"%len(atoms))
  for n in range(len(atoms)):
    fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.close()
  return True

def TXYZ2XYZi(TXYZ):
  data = readTXYZ(TXYZ)
  atoms = data[0]; coord = data[1]
  fileName = TXYZ.split(".txyz")[0]
  fout=open(fileName+".xyz",'w')
  fout.write("%3s\n#Generated by chemFileConvert.py\n"%(len(atoms)-1))
  for n in range(len(atoms)-1):
    fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.close()
  return True

def TXYZ2TXYZ(TXYZ):
  data = readTXYZ(TXYZ)
  atoms = data[0];coord = data[1]
  order = data[2];types = data[3]
  connections = data[4]
  number = int(order[-1])-len(order)
  for i in range(len(order)):
    order[i]=int(order[i])-number

  for i in range(len(connections)):
    for j in range(len(connections[i])):
      connections[i][j]=int(connections[i][j])-number

  fileName = TXYZ.split(".xyz")[0]
  fout=open(fileName+".xyz_1",'w')
  fout.write("%3s #Generated by chemFileConvert.py\n"%len(atoms))
  for n in range(len(atoms)):
    if len(connections[n])==0:
      fout.write("%6s%3s%14.7f%14.7f%14.7f%6s\n"%(order[n],atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]), types[n]))
    if len(connections[n])==1:
      fout.write("%6s%3s%14.7f%14.7f%14.7f%6s%3s\n"%(order[n],atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]), types[n], connections[n][0]))
    if len(connections[n])==2:
      fout.write("%6s%3s%14.7f%14.7f%14.7f%6s%3s%3s\n"%(order[n],atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]), types[n], connections[n][0], connections[n][1]))
    if len(connections[n])==3:
      fout.write("%6s%3s%14.7f%14.7f%14.7f%6s%3s%3s%3s\n"%(order[n],atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]), types[n], connections[n][0], connections[n][1], connections[n][2]))
    if len(connections[n])==4:
      fout.write("%6s%3s%14.7f%14.7f%14.7f%6s%3s%3s%3s%3s\n"%(order[n],atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]), types[n], connections[n][0], connections[n][1], connections[n][2], connections[n][3]))
  fout.close()
  #os.system("mv %s.xyz %s"%(fileName, TXYZ))
  return True

def COM2TXYZ(COM,conectFile):
  data = readCOM(COM)
  atoms = data[0]
  coord = data[1]
  oFile = open("tmp","w")
  oFile.write("%3s #Generated by chemFileConvert.py\n"%len(atoms))
  for i in range(len(atoms)):
    oFile.write("%3s%3s%14.7f%14.7f%14.7f\n"%(str(i+1), atoms[i], coord[i][0], coord[i][1], coord[i][2]))
  oFile.close()
  os.system("paste tmp %s > %s"%(conectFile, COM.split(".")[0]+".txyz")) 
  return True

def TXYZ2TAIL(TXYZ):
  data = readTXYZ(TXYZ)
  atoms = data[0];coord = data[1]
  order = data[2];types = data[3]
  connections = data[4]
  number = int(order[-1])-len(order)
  for i in range(len(order)):
    order[i]=int(order[i])-number

  for i in range(len(connections)):
    for j in range(len(connections[i])):
      connections[i][j]=int(connections[i][j])-number

  fileName = TXYZ.split(".xyz")[0]
  fout=open(fileName+".tail",'w')
  fout.write("     \n")
  for n in range(len(atoms)):
    if len(connections[n])==0:
      fout.write("%6s\n"%(types[n]))
    if len(connections[n])==1:
      fout.write("%6s%3s\n"%(types[n], connections[n][0]))
    if len(connections[n])==2:
      fout.write("%6s%3s%3s\n"%(types[n], connections[n][0], connections[n][1]))
    if len(connections[n])==3:
      fout.write("%6s%3s%3s%3s\n"%(types[n], connections[n][0], connections[n][1], connections[n][2]))
    if len(connections[n])==4:
      fout.write("%6s%3s%3s%3s%3s\n"%(types[n], connections[n][0], connections[n][1], connections[n][2], connections[n][3]))
  fout.close()
  return True

def XYZ2QCHEM(XYZ):
  atoms = readXYZ(XYZ)[0]
  coord = readXYZ(XYZ)[1]
  
  qfile = open(XYZ.split(".xyz")[0]+".qchem", "w")
  qfile.write("$molecule\n")
  qfile.write("%s\n"%chrgMult)
  for i in range(len(atoms)):
    qfile.write("%3s%12.6f%12.6f%12.6f\n"%(atoms[i], coord[i][0], coord[i][1], coord[i][2]))
  qfile.write("$end\n\n")
  qfile.write("$rem\n")
  qfile.write(" Basis = %s\n"%basis)
  qfile.write(" JOB_TYPE = %s\n"%jobType)
  qfile.write(" METHOD = %s\n"%QM)
  qfile.write("$end\n\n")
  qfile.close()
    

def ARC2ARC(ARC):
  boxInfo = "90.000000   90.000000   90.000000"
  fileName = ARC.split(".arc")[0]
  fout=open(fileName+".arc_1",'w')
  lines = file(ARC).readlines()
  for line in lines:
    if boxInfo not in line:
      fout.write(line)
  fout.close()
  return True

def XYZ2TXYZ(dimerXYZ, mono1TXYZ, mono2TXYZ):
  dimerCoord = []
  mono1 = []
  mono2 = []
  #standard xyz file, with second line as comment
  for line in open(dimerXYZ).readlines()[2:]:
    dimerCoord.append(line.split("\n")[0])
  #txyz file, no comment in the second line
  for line in open(mono1TXYZ).readlines()[1:]:
    data = line.split()
    mono1.append([data[0],] + data[5:])
  #txyz file, no comment in the second line
  for line in open(mono2TXYZ).readlines()[1:]:
    data = line.split()
    mono2.append([data[0],] + data[5:])
  natom = len(dimerCoord)
  filename = dimerXYZ.split(".xyz")[0]+".txyz"
  ofile = open(filename, "w")
  ofile.write("%3s Generated by chemFileConvert.py\n" %natom)
  for i in range(len(mono1)):
    if len(mono1[i]) == 2:
      ofile.write( "%3s %s %3s\n" %(mono1[i][0], dimerCoord[i], mono1[i][1]))
    elif len(mono1[i]) == 3:
      ofile.write( "%3s %s %3s%3s\n" %(mono1[i][0], dimerCoord[i], mono1[i][1], mono1[i][2]))
    elif len(mono1[i]) == 4:
      ofile.write( "%3s %s %3s%3s%3s\n" %(mono1[i][0], dimerCoord[i], mono1[i][1], mono1[i][2], mono1[i][3]))
    elif len(mono1[i]) == 5:
      ofile.write( "%3s %s %3s%3s%3s%3s\n" %(mono1[i][0], dimerCoord[i], mono1[i][1], mono1[i][2], mono1[i][3],mono1[i][4]))
    elif len(mono1[i]) == 6:
      ofile.write( "%3s %s %3s%3s%3s%3s%3s\n"%(mono1[i][0], dimerCoord[i], mono1[i][1], mono1[i][2], mono1[i][3],mono1[i][4], mono1[i][5]))
    elif len(mono1[i]) == 7:
      ofile.write( "%3s %s %3s%3s%3s%3s%3s\n"%(mono1[i][0], dimerCoord[i], mono1[i][1], mono1[i][2], mono1[i][3],mono1[i][4], mono1[i][5],mono1[i][6]))
  for i in range(len(mono2)):
    if len(mono2[i]) == 2:
      ofile.write( "%3s %s %3s\n" %(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1]))
    elif len(mono2[i]) == 3:
      ofile.write( "%3s %s %3s%3s\n" %(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1], len(mono1) + int(mono2[i][2])))
    elif len(mono2[i]) == 4:
      ofile.write( "%3s %s %3s%3s%3s\n" %(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1], len(mono1) + int(mono2[i][2]), len(mono1) + int(mono2[i][3])))
    elif len(mono2[i]) == 5:
      ofile.write( "%3s %s %3s%3s%3s%3s\n" %(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1], len(mono1) + int(mono2[i][2]), len(mono1) + int(mono2[i][3]),len(mono1) + int(mono2[i][4])))
    elif len(mono2[i]) == 6:
      ofile.write( "%3s %s %3s%3s%3s%3s%3s\n"%(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1], len(mono1) + int(mono2[i][2]), len(mono1) + int(mono2[i][3]),len(mono1) + int(mono2[i][4]), len(mono1) + int(mono2[i][5])))
    elif len(mono2[i]) == 7:
      ofile.write( "%3s %s %3s%3s%3s%3s%3s\n"%(len(mono1) + int(mono2[i][0]), dimerCoord[i], mono2[i][1], len(mono1) + int(mono2[i][2]), len(mono1) + int(mono2[i][3]),mono2[i][4], len(mono1) + int(mono2[i][5],mono2[i][6])))
  return 
