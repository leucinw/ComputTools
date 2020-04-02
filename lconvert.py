
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


#  TXYZ2COM:  From TINKER xyz file to Gaussian input file
#   LOG2COM:   From Gaussian output file to Gaussian input file
#   COM2COM:   From Gaussian input file to new input file
#   XYZ2COM:   From common xyz file to Gaussian input file

#   COM2TXYZ:  From Gaussian input to TINKER XYZ file
#  TXYZ2TXYZ: From TINKER xyz file to a txyz file that all atoms are in order
#    ARC2ARC:   From TINKER traj file to traj file
#   XYZ2TXYZ:  From common XYZ to Tinker XYZ file

#  TXYZ2XYZ:   From TINKER xyz file to common XYZ file (most software readable)
#   COM2XYZ:   From Gaussian input file to common XYZ file
#   LOG2XYZ:   From Gaussian output file to normal xyz file
#   XYZ2XYZ:   From normal xyz file to normal xyz file in different distances
#  XYZ2MONO:   From dimer xyz file to monomer xyz file

#  XYZ2QCHEM: From common XYZ file to Qchem file
#  XYZ2PSI4:  From normal xyz file to psi4 file(SAPT calculation)
#  XYZ2CBS:   From normal xyz file to psi4 file(CBS calculation)
#  XYZ2ORCA:  From normal xyz file to ORCA input file 
#  XYZ2CON:   From normal xyz file extract conection information
#  SDF2PNG:   From sdf file to 2d structure image

import os
from readChemFile import *
from operation import *
from chem import *

# 1. com to com
def COM2COM(COM):
  atoms, coord = readCOM(COM)
  fileName = COM.replace(".com", "com_1")
  with open(fileName, 'w') as fout: 
    fout.write("%chk="+"%s.chk\n"%fileName.split(".com")[0])
    fout.write("%nproc="+nproc+"\n%Mem=32GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
    for n in range(len(atoms)):
      fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
    fout.write("\n")
  return

# 2. common xyz to com
def XYZ2COM(XYZ, method="HF", basis="3-21G ", nproc="8", charge="0", multiplicity="1", jobtype="SP", memory="100GB", extra= " ", restraint = False):
  chrgMult=charge + " " + multiplicity
  jobType = jobtype 
  extraKeyword = " MaxDisk=100GB IOP(5/13=1)" + extra 
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  with open(fileName + ".com", 'w') as fout:
    fout.write("%chk="+"%s.chk\n"%fileName)
    fout.write("%nproc="+nproc+"\n%Mem="+memory+"\n#p "+method+"/"+basis + extraKeyword + " %s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
    for n in range(len(atoms)):
      fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
    fout.write("\n")
    if restraint==True:
      fout.write("* 1 2 * R\n")
      fout.write("* 2 11 * R\n")
      fout.write("4 1 2 11 F\n")
      fout.write("1 2 11 12 F\n")
  return 

def XYZ2COM_BSSE(XYZ, nfirst, chrgMult, optCartesian=False):
  #QM='CCSD(T)'
  QM='mp2'
  basis ='genECP '
  nproc     ='8'
  jobType = "SP "
  if optCartesian == True:
    jobType = "Opt=(tight,Cartesian) "
  extraKeyword = " counterpoise=2 IOP(5/13=1) Symmetry=None MaxDisk=100GB 5d 7f Pseudo=Read "
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  fout=open(fileName+".com",'w')
  fout.write("%NoSave\n")
  fout.write("%chk="+"%s.chk\n"%fileName)
  fout.write("%nproc="+nproc+"\n%Mem=64GB\n#p "+QM+"/"+basis+extraKeyword+"%s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(nfirst):
    fout.write("%3s(Fragment=1)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  for n in range(nfirst,len(atoms)):
    fout.write("%3s(Fragment=2)             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

def XYZ2PSI4(XYZ, nfirst, chrgMult_1, chrgMult_2, basis, saptDFT=False):
  data = readXYZ(XYZ)
  atoms = data[0]
  coord = data[1]
  fout=open(XYZ.split(".xyz")[0]+".psi4","w")
  fout.write("#psi4\n#\nmemory 48 gb\n\nmolecule  {\n%s\n"%chrgMult_1)
  for n in range(nfirst):
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("--\n%s\n"%chrgMult_2)
  for n in range(nfirst, len(coord)):
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("units angstrom\nno_reorient  #important for SAPT in psi4, default\n")
  fout.write("symmetry c1  #important for SAPT in psi4, default\n}\n\n")
  if basis == "TZ":
    fout.write("set {\nbasis aug-cc-pVTZ\nscf_type DF\n")
    if saptDFT==True:
      fout.write("DO_IND_EXCH_SINF True\nSAPT_DFT_FUNCTIONAL HF\n}\nenergy('sapt(dft)')\n")
    else:
      fout.write("freeze_core True\n}\nenergy('sapt2+')\n")
  if basis == "DZ":
    fout.write("set {\nbasis aug-cc-pVDZ\nscf_type DF\n")
    if saptDFT==True:
      fout.write("DO_IND_EXCH_SINF True\nSAPT_DFT_FUNCTIONAL HF\n}\nenergy('sapt(dft)')\n")
    else:
      fout.write("freeze_core True\n}\nenergy('sapt2+')\n")

  return True

def XYZ2CBS(XYZ, chrgMult = "0 1"):
  atoms, coord = readXYZ(XYZ)
  fout=open(XYZ.split(".xyz")[0]+".psi4","w")
  fout.write("#psi4\n#\nmemory 20 gb\n\nmolecule  {\n%s\n"%chrgMult)
  for n in range(len(atoms)):
    fout.write("%3s   %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("units angstrom\nno_reorient  #important for SAPT in psi4, default\n")
  fout.write("symmetry c1  #important for SAPT in psi4, default\n}\n\n")
  fout.write("set {\nscf_type DF\nmp2_type DF\ne_convergence 7\nreference rhf\n}\n")
  fout.write("energy('mp2/aug-cc-pv[dt]z')\n")
  #fout.write("energy(cbs, corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')\n")
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
  data = readLOG(LOG)
  nproc = "8"
  QM = "MP2"
  basis = "6-31G*"
  extraKeyword = " "
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

def XYZ2TXYZ(XYZ,indexInit):
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  connections = XYZ2CON(XYZ)
  fileName = XYZ.split(".xyz")[0]
  natom = len(atoms)
  ofile = open(fileName+".txyz", "w")
  ofile.write("%3s Generated by chemFileConvert.py\n" %natom)
  print(len(atoms))
  print(natom)
  print(len(coord))
  print(len(connections))
  for i in range(natom):
    ofile.write("%3s %3s %12.6f%12.6f%12.6f %3s %s\n" %(i+1, atoms[i], coord[i][0], coord[i][1], coord[i][2], (int(indexInit) + i), connections[i]))
  return

def XYZ2TXYZ_Water(XYZ):
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  natom = len(atoms)
  ofile = open(fileName+".txyz", "w")
  ofile.write("%3s Generated by chemFileConvert.py\n" %natom)
  for i in range(natom):
    if i%3==0: #Oxygen
      ofile.write( "%3s %3s %12.6f%12.6f%12.6f 1 %3s %3s\n" %(i+1, atoms[i], coord[i][0], coord[i][1],coord[i][2], i+2, i+3))
    if i%3==1: #Hydrogen 1
      ofile.write( "%3s %3s %12.6f%12.6f%12.6f 2 %3s\n" %(i+1, atoms[i], coord[i][0], coord[i][1],coord[i][2], i))
    if i%3==2: #Hydrogen 2 
      ofile.write( "%3s %3s %12.6f%12.6f%12.6f 2 %3s\n" %(i+1, atoms[i], coord[i][0], coord[i][1],coord[i][2], i-1))
  return

def XYZ2ORCA(XYZ):
  data = readXYZ(XYZ)
  atoms = data[0];coord = data[1]
  fileName = XYZ.split(".xyz")[0]
  fout=open(fileName+".com",'w')
  fout.write("%chk="+"%s.chk\n"%fileName)
  fout.write("%nproc="+nproc+"\n%Mem="+memory+"\n#p "+method+"/"+basis+" %s\n\n%s \n\n%s\n"%(jobType,jobType,chrgMult))
  for n in range(len(atoms)):
    fout.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))
  fout.write("\n")
  fout.close()
  return True

def XYZ2CON(xyz):
  connection = []
  os.system("babel -ixyz %s -opdb tmp.pdb"%xyz)
  lines = readwholefile("tmp.pdb")
  for line in lines:
    if "conect" in line:
      connection.append(line[11:-1])
  return connection

def SDF2PNG(SDF):
  # make sure you can run obabel
  cmdstr = "obabel -isdf %s -O %s.png -p 300"%(SDF, SDF.split(".")[0])
  os.system(cmdstr)
  return True
