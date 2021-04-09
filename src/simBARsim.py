

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


# Integrate BAR simulations in one script

import os
import sys
import time
import numpy as np
import subprocess

# color
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\33[93m'
ENDC = '\033[0m'

def sub(node, runfile, userpath=None):
  if userpath:
    cwd = userpath
  else: 
    cwd = os.getcwd()
  cmdstr = 'ssh %s "cd %s; nohup sh %s 2>err.log &" &' % (node, cwd, runfile)
  subprocess.run(cmdstr, shell=True)
  return

def readini():
  argDict = {}
  lines = [line for line in open("bar.ini").readlines()]
  for line in lines:
    if ("#" not in line) and (line != "\n"):
      setting = line.split("\n")[0].split("=")
      key = setting[0].upper()
      val = setting[1]
      argDict[key] = val
  return argDict

def main():
  argDict = readini()
  orderparams = [] 
  nodelist = []

  if argDict["ORDERPARAMS"] == "read":
    for line in open("orderparams").readlines():
      if "#" not in line:
        orderparams.append(line.split())
  else:
    for prm in argDict["ORDERPARAMS"].split(">"):
      orderparams.append(prm.split(","))

  if argDict["NODELIST"]=="read":
    lines = open("nodelist").readlines()
    for line in lines:
      if "#" not in line[0]:
        nodelist.append(line.split("\n")[0])
  else:
    for node in argDict["NODELIST"].split(","):
      nodelist.append(node)

  fname = argDict["DIRNAME"]
  currDir = os.getcwd()

  def setup():
    for i in range(len(orderparams)):
      elb = orderparams[i][0]
      vlb = orderparams[i][1]
      dirname = "%s-%s-%s"%(argDict["DIRNAME"], "%03d"%int(float(elb)*100), "%03d"%int(float(vlb)*100))      
      if not os.path.isdir(dirname): 
        os.mkdir(dirname)
      subprocess.run("cp %s temp.key"%argDict["TINKERKEY"], shell=True)
      with open("temp.key","a") as temp:
        temp.write("ele-lambda %10.3f\n"%float(elb))
        temp.write("vdw-lambda %10.3f\n"%float(vlb))
      copystr = "cp %s %s %s"%(argDict["TINKERXYZ"], argDict["DYNAMICRUN"], dirname)
      subprocess.run(copystr, shell=True)
      subprocess.run("mv temp.key %s/%s"%(dirname, argDict["TINKERKEY"]), shell=True)
      if i < len(orderparams)-1:  
        elb_ = orderparams[i+1][0]
        vlb_ = orderparams[i+1][1]
        dirname_ = "%s-%s-%s"%(argDict["DIRNAME"], "%03d"%int(float(elb_)*100), "%03d"%int(float(vlb_)*100))      
        lines = open(argDict["BARRUN"]).readlines()
        for j in range(len(lines)):
          if "export nextstep" in lines[j]:
            lines[j] = "export nextstep=%s\n"%dirname_
        with open("temp.sh", "w") as temp:
          for line in lines:
            temp.write(line)
        subprocess.run("mv temp.sh %s/%s"%(dirname, argDict["BARRUN"]),shell=True)
    print(GREEN + "BAR simulation folders generated!" + ENDC)
    return
  
  def dynamic():
    """submit the dynamic job of each window"""
    for i in range(len(orderparams)):
      tmp = orderparams[i]
      dirname = os.path.join(currDir, fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
      if os.path.isfile(os.path.join(dirname, argDict["TINKERXYZ"].replace(".xyz",".arc"))):
        print(YELLOW + " .arc file exists in %s !"%dirname + ENDC)
      else:
        sub(nodelist[i], argDict["DYNAMICRUN"], dirname) 
        print(GREEN + " Submitted dynamic job on %s"%nodelist[i] + ENDC)

    """check the completeness of dynamic job"""
    while True: 
      nfinish = 0
      for i in range(len(orderparams)):
        tmp = orderparams[i]
        dirname = os.path.join(currDir, fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
        f = 'liquid.log'
        f = os.path.join(dirname,f)
        if os.path.isfile(f):
          tail = subprocess.check_output("tail -n10 %s"%f, shell=True).decode("utf-8")
          line = ''.join(list(tail))
          if "ns/day" in line:
            nfinish += 1
      if nfinish == len(orderparams):
        print(GREEN + " All DYNAMIC jobs have finished, cheers!" + ENDC)
        break
      else:
        print(RED + " Waiting for DYNAMIC jobs to finish!" + ENDC)
        time.sleep(30.0)
    return
  
  def bar():
    """submit bar simulations after trajectories generated"""
    for i in range(len(orderparams)-1):
      tmp = orderparams[i]
      dirname = os.path.join(currDir, fname+"-%03d-%03d"% ((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
      if os.path.isfile(os.path.join(dirname, argDict["TINKERXYZ"].replace(".xyz", ".bar"))):
        subprocess.run("rm -f %s/freeEnergy.log %s"%(dirname,os.path.join(dirname, argDict["TINKERXYZ"].replace(".xyz",".bar"))),shell=True)
        sub(nodelist[i], argDict["BARRUN"], dirname)
        print(YELLOW + " Submitted BAR job on %s"%nodelist[i] + ENDC)
      else:
        sub(nodelist[i], argDict["BARRUN"], dirname) 
        print(GREEN + " Submitted BAR job on %s"%nodelist[i] + ENDC)

    """check the completeness of bar job"""
    while True: 
      nfinish = 0
      for i in range(len(orderparams)-1):
        tmp = orderparams[i]
        dirname = os.path.join(currDir, fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
        f = 'freeEnergy.log'
        f = os.path.join(dirname,f)
        if os.path.isfile(f):
          tail = subprocess.check_output("tail -n100 %s"%f, shell=True).decode("utf-8")
          line = ''.join(list(tail))
          if "Free Energy via BAR Bootstrap" in line:
            nfinish += 1
      if nfinish == (len(orderparams)-1):
        print(GREEN + " All BAR jobs have finished, cheers!" + ENDC)
        break
      else:
        print(RED + " Waiting for BAR jobs to finish !" + ENDC)
        time.sleep(30.0)
    return
  
  def result():
    Free = []
    Err = []
    Dir = []
    fo = open("result.txt", "w")
    for i in range(len(orderparams)-1):
      tmp = orderparams[i]
      dirname = os.path.join(currDir, fname+"-%03d-%03d"% ((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
      Dir.append(fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
      if os.path.isfile(os.path.join(dirname, "freeEnergy.log")):
        for line in open(os.path.join(dirname,"freeEnergy.log")).readlines():
          if "Free Energy via BAR Iteration" in line:
            Free.append(float(line.split()[-4]))
            Err.append(float(line.split()[-2]))
      else:
            print(RED + "Free Energy calculation not complete in %s"%fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))) + ENDC)
    [l10,l20] = orderparams[0]
    [l11,l21] = orderparams[1]
    if (float(l10)+float(l20)) > (float(l11)+float(l21)):
      disappear = True
    else:
      disappear = False
    if disappear:
      for i,j,k in zip(Dir,Free,Err):
        print(YELLOW + "%20s%15.8f%15.8f"%(i,-j,k) + ENDC)
        fo.write("%20s%15.8f%15.8f\n"%(i,-j,k))
      totFree = np.array(Free).sum()
      totErr = np.sqrt(np.square(np.array(Err)).sum())
      print(GREEN + "%20s%15.8f%15.8f"%("Total Free Energy", -totFree, totErr) + ENDC)
      fo.write("%20s%15.8f%15.8f\n"%("Total Free Energy", -totFree, totErr))
    else:
      for i,j,k in zip(Dir,Free,Err):
        print(YELLOW + "%20s%15.8f%15.8f"%(i,j,k) + ENDC)
        fo.write("%20s%15.8f%15.8f\n"%(i,j,k))
      totFree = np.array(Free).sum()
      totErr = np.sqrt(np.square(np.array(Err)).sum())
      print(GREEN + "%20s%15.8f%15.8f"%("Total Free Energy", totFree, totErr) + ENDC)
      fo.write("%20s%15.8f%15.8f\n"%("Total Free Energy", totFree, totErr))
    return

  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result}
  for action in sys.argv[1:]:
    if action in actions:
      actions[action]()
    else:
      sys.exit(RED + "Usage: python simBARsim.py setup/dynamic/bar/result" + ENDC)

if len(sys.argv) == 1:
  sys.exit(RED + "Usage: python simBARsim.py setup/dynamic/bar/result" + ENDC)
else:
  main()
