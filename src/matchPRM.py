

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Generate tinker prm file for new atomtype txyz

''' python matchPRM.py template.prm template.txyz dealwith.txyz energyterm'''

import os,sys,subprocess

tinkerexedir = "$TINKER8C8"

def readTXYZ(TXYZ, singleAtom=False):
  atoms=[];coord=[]; tails = []
  order=[];types=[];connections=[]
  for line in open(TXYZ).readlines()[1:]: 
    data=line.split()
    order.append(data[0])
    types.append(data[5])
    connections.append(data[6:])
    if singleAtom:
      atoms.append(data[1][0])
    else:
      atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
  return atoms,coord,order,types,connections

def matchAllTerms(templateprm, templatexyz, dealwithxyz, term):
  #
  # define some dictionary for later use 
  #
  termDict = {"bond":"BOND.prm", "angle":"ANGLE.prm", "torsion":"TORSION.prm", \
              "strbnd":"STRBND.prm", "opbend":"OPBEND.prm", "pitors":"PITORS.prm"}
  termKeyDict = {"bond":"BONDterm", "angle":"ANGLEterm", "torsion":"TORSIONterm", \
                 "strbnd":"STRBNDterm", "opbend":"OPBENDterm", "pitors":"PITORSterm"}
  termOutDict = {"bond": "Bond Stretching Parameters :", \
                 "angle": "Angle Bending Parameters :", \
                 "torsion": "Torsional Angle Parameters :", \
                 "strbnd": "Stretch-Bend Parameters :",\
                 "pitors": "Pi-Orbital Torsion Parameters :",\
                 "opbend": "Out-of-Plane Bend Parameters :"}
  #
  # for 3-fold torsion
  #
  onefold   = "0.000    0.0   1   "
  twofold   = "0.000  180.0   2   "
  threefold = "0.000    0.0   3   "

  #
  # make an atom-to-type dictionary
  #
  _,_,_,types,_ = readTXYZ(dealwithxyz)
  typesDict = {}
  for i in range(len(types)):
    typesDict[i+1] = types[i]

  #
  # call tinker analyze.x to get parameters
  #
  with open("tinker.key",'w') as f:
    f.write("parameters %s\n"%templateprm)
    f.write("%s only\n"%termKeyDict[term.lower()])
  cmdstr = "%s/analyze.x %s -k tinker.key P > tmp.out "%(tinkerexedir, templatexyz)   
  subprocess.run(cmdstr,shell=True)

  #
  # find the parameters
  #
  lines = open("tmp.out").readlines()
  for line in lines:
    if termOutDict[term.lower()] in line:
      idx = lines.index(line)
  idx +=4

  #
  # map the parameters back to new types
  #
  repeated = []
  ttt = open(termDict[term.lower()], "w")
  for i in range(idx, len(lines), 1):
    if term.lower() == "bond":
      a1, a2, ks, bnd = lines[i].split()[1:5]
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      prmstr = ' '.join([ks, bnd])
      atomliststr = 2*"%5s"%(a1_, a2_)
    elif term.lower() == "pitors":
      a1, a2, ks = lines[i].split()[1:4]
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      prmstr = ks 
      atomliststr = 2*"%5s"%(a1_, a2_)
    elif term.lower() == "angle":
      a1, a2, a3, kb, ang = lines[i].split()[1:6]
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      a3_ = typesDict[int(a3)]
      prmstr = ' '.join([kb, ang])
      atomliststr = 3*"%5s"%(a1_, a2_, a3_)
    elif term.lower() == "strbnd":
      a1, a2, a3, ksb1, ksb2 = lines[i].split()[1:6]
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      a3_ = typesDict[int(a3)]
      prmstr = ' '.join([ksb1, ksb2])
      atomliststr = 3*"%5s"%(a1_, a2_, a3_)
    elif term.lower() == "opbend":
      a1, a2, a3, a4, kopb = lines[i].split()[1:6]
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      a3_ = typesDict[int(a3)]
      a4_ = typesDict[int(a4)]
      atomliststr = 4*"%5s"%(a1_, a2_, a3_, a4_)
      prmstr = ' '.join([kopb])
    elif term.lower() == "torsion":
      a1, a2, a3, a4 = lines[i].split()[1:5]
      if "0/1" in lines[i]:
        index = int(lines[i].split().index("0/1")) - 1
        onefold_ = '   '.join([lines[i].split()[index], '0.0', '1'])
      else:
        onefold_ = onefold
      if "180/2" in lines[i]:
        index = int(lines[i].split().index("180/2")) - 1
        twofold_ = '   '.join([lines[i].split()[index], '180.0', '2'])
      else:
        twofold_ = twofold
      if "0/3" in lines[i]:
        index = int(lines[i].split().index("0/3")) - 1
        threefold_ = '   '.join([lines[i].split()[index], '0.0', '3'])
      else:
        threefold_ = threefold
      a1_ = typesDict[int(a1)]
      a2_ = typesDict[int(a2)]
      a3_ = typesDict[int(a3)]
      a4_ = typesDict[int(a4)]
      atomliststr = 4*"%5s"%(a1_, a2_, a3_, a4_)
      prmstr = ' '.join([onefold_, twofold_, threefold_])
    else:
      print("%s currently not supported!")
      break
    #
    # write out parameters    
    #
    writingstr = "%s %s %s\n"%(term, atomliststr, prmstr)
    if writingstr not in repeated:
      repeated.append(writingstr)
      ttt.write(writingstr)
  ttt.close()       
  return

def main():
  templateprm = sys.argv[1]
  templatexyz = sys.argv[2]
  dealwithxyz = sys.argv[3]
  energyterm  = sys.argv[4]
  matchAllTerms(templateprm, templatexyz, dealwithxyz, energyterm)
  return

if __name__ == "__main__":
  main()
