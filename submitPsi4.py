
#===================================
# 	     Chengwen Liu              #
# 	   liuchw2010@gmail.com        #
# 	University of Texas at Austin  #
#===================================


import os
import sys
import time
from readDataFile import *

''' 
    Check the availability of a certain node 
    Return True/False
'''
def checkPsi4Jobs(node):
  psi4Jobs = 0
  cmdstr = "ssh %s 'top -n1 -b' >top.tmp" % node
  os.system(cmdstr)
  cmdstr = "ssh %s 'ps aux | grep liuchw | grep psi4 ' >ps.tmp" % node
  os.system(cmdstr)
  psResult = readWholeFile("ps.tmp")
  topResult = readWholeFile("top.tmp")
  for ps in psResult:
    if "psi4" in ps:
      psi4Jobs += 1
  for top in topResult:
    if "psi4" in top:
      psi4Jobs += 1

  if psi4Jobs == 0:
    print("%s is IDLE, ready to submit your jobs!"%node)
    return True 
  else:
    print("%s is BUSY, come back later!"%(node))
    return False

'''
    Submit one psi4 job on an idle node 
'''
def sub_one_sapt(node, psi4):
  fname = psi4.split(".psi4")[0]
  cwd = os.getcwd()
  cmdstr = 'ssh %s "source ~/.bashrc.psi4conda; cd %s; nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.log &" &' % (node, cwd, fname, fname)
  os.system(cmdstr)
  return

''' auto-submit jobs on specified nodelist
    for specified files in filelist 
    1. check idle nodes and save to an idlelist
    2. submit jobs on idle nodes
    3. keep doing 1,2 until all jobs are submitted 
'''
def sub_multiple_jobs(nodelist, filelist):
  currentJob = 0
  nodes = readWholeFile(nodelist)
  files = readWholeFile(filelist)
  while(True):
    idleNodes = []
    if currentJob == (len(files)):
      break
    for eachNode in nodes:
      if checkPsi4Jobs(eachNode) == True:
        idleNodes.append(eachNode)
    if idleNodes == []:
      print("So I will sleep for 10 minites!")
      time.sleep(60.0*10.0)
    else: 
      for eachNode in idleNodes:
        if currentJob == (len(files)):
          break
        sub_one_sapt(eachNode, files[currentJob])
        print("So I submitted %s on %s!"%(files[currentJob], eachNode))
        currentJob += 1
  return 

sub_multiple_jobs("nodelist", "filelist")
