

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# After poltype run, we have ttt.xyz and ttt.key as the final force field txyz and parameter file
# It is very often that files have their class definition overlapped.
# Here I will use this script to automatically deal with this 
# This should be run at the top level of a bunch of completed poltype jobs
#

# 1. walk through the directory and find ttt.xyz and ttt.key
# 2. apply the new atom classes to ttt.xyz and ttt.key  
# 3. output xyz files and appended key file to a new diredctory

import os,sys
from readDataFile import *
from readChemFile import *
from writeChemFile import *

def modOneTXYZ(TXYZ, classNum, outXYZ):
  data = readTXYZ(TXYZ)
  atoms = data[0];coord = data[1]
  order = data[2];types = data[3]
  connections = data[4]
  types = [int(types[i]) for i in range(len(types))]
  minType = min(types)
  for i in range(len(types)):
    types[i] = (types[i] - minType + classNum)
  writeTXYZ(outXYZ, order, atoms, coord, types, connections)
  return len(set(types))

def modOnePRM(PRM, classNum, outPRM, term):
  # atom
  if term == "atom ":
    data = readTinkerPRM(PRM, "atom ")
    atoms = data[0]
    atomList = data[1]
    typeList = data[2]
    remainList = data[3]
    typeList = [int(typeList[i]) for i in range(len(typeList))]
    minType = min(typeList)
    for i in range(len(typeList)):
      typeList[i] = (typeList[i] - minType + classNum)
      atomList[i] = (atomList[i] - minType + classNum)
    ofile = open(outPRM, "a")  
    for i in range(0, len(atoms)):
      remainStr = ''
      for remain in remainList[i]:
        remainStr += (remain + "   ")
      ofile.write("%4s %3s %3s %s\n"%(atoms[i], atomList[i], typeList[i], remainStr))
    ofile.close()
    return len(set(typeList))

  # polarize
  if term == "polarize ":
    data = readTinkerPRM(PRM, "polarize ")
    polarizes = data[0]
    atomList = data[1]
    atomicPolarList = data[2]
    dampingList = data[3]
    polarGroupList = data[4]

    atomList = [int(atomList[i]) for i in range(len(atomList))]
    minType = min(atomList)
    for i in range(len(atomList)):
      atomList[i] = (atomList[i] - minType + classNum)
      for j in range(len(polarGroupList[i])):
        polarGroupList[i][j] = (polarGroupList[i][j] - minType + classNum)

    ofile = open(outPRM, "a")  
    for i in range(0, len(polarizes)):
      polarGrpStr = ''
      for polGrp in polarGroupList[i]:
        polarGrpStr += (str(polGrp) + " ")
      ofile.write("%8s %5s %8s %8s %s\n"%(polarizes[i], atomList[i], atomicPolarList[i], dampingList[i], polarGrpStr))
    ofile.close()
    return len(set(atomList))

#os.system("rm -rf tinker.key")
#classNumber = 1 
#for eachFile in readWholeFile("filelist"):
#  addedNum = modOnePRM(eachFile+"/ttt.key", classNumber, "tinker.key", "atom ")
#  classNumber += addedNum
#
#classNumber = 1 
#for eachFile in readWholeFile("filelist"):
#  addedNum = modOnePRM(eachFile+"/ttt.key", classNumber, "tinker.key", "polarize ")
#  classNumber += addedNum
#  print(classNumber)

classNumber = 1
os.system("rm -rf myXYZs/*.xyz")
for eachFile in readWholeFile("filelist"):
  print(os.listdir(eachFile))
  if "ttt.xyz" in os.listdir(eachFile):
    filename = eachFile+"/ttt.xyz"
  else:  
    filename = eachFile+"/ttt.txyz"
  addedNum = modOneTXYZ(filename, classNumber, "temp.xyz")
  os.system("mv temp.xyz ./myXYZs/%s.txyz"%eachFile)
  classNumber += addedNum

