

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


# Integrate BAR simulations in one script

import os,sys
import numpy as np
import subprocess

def sub(node, input_file, userpath=None):
  fname, ext = os.path.splitext(input_file)
  if userpath:
    cwd = userpath
  else: 
    cwd = os.getcwd()
  cmdstr = 'ssh %s "cd %s; nohup sh %s.sh 2>err.log &" &' % (node, cwd, fname)
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

def setup():
  argDict = readini()
  orderparams = [] 
  if argDict["ORDERPARAMS"] == "read":
    lines = open("orderparams").readlines()
    for line in lines:
      if "#" not in line:
        orderparams.append(line.split())
  else:
    for prm in argDict["ORDERPARAMS"].split(">"):
      orderparams.append(prm.split(","))
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
    copystr = "cp %s %s %s %s"%(argDict["TINKERXYZ"], argDict["TINKERPRM"], argDict["DYNAMICRUN"], dirname)
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
  return

def dynamic():
  currDir = os.getcwd()
  argDict = readini()
  orderparams = []
  nodelist = []
  orderparams = [] 
  fname = argDict["DIRNAME"]
  if argDict["ORDERPARAMS"] == "read":
    lines = open("orderparams").readlines()
    for line in lines:
      if "#" not in line:
        orderparams.append(line.split())
  else:
    for prm in argDict["ORDERPARAMS"].split(">"):
      orderparams.append(prm.split(","))
  if argDict["NODELIST"]=="read":
    nodelist = [line.split("\n")[0] for line in open("nodelist").readlines()]
  else:
    for node in argDict["NODELIST"].split(","):
      nodelist.append(node)
  for i in range(len(orderparams)):
    tmp = orderparams[i]
    dirname = currDir + "/" + fname + "-%03d-%03d"% ((int(float(tmp[0])*100)), (int(float(tmp[1])*100)))
    if os.path.isfile(os.path.join(dirname, argDict["TINKERXYZ"].split(".xyz")[0]+".arc")):
      print(" .arc file exists in %s !"%dirname)
    else:
      sub(nodelist[i], argDict["DYNAMICRUN"], dirname) 
      print(" Submitted dynamic job on %s !"%nodelist[i])
  return

def bar():
  currDir = os.getcwd()
  argDict = readini()
  orderparams = []
  nodelist = []
  orderparams = [] 
  fname = argDict["DIRNAME"]
  if argDict["ORDERPARAMS"] == "read":
    lines = open("orderparams").readlines()
    for line in lines:
      if "#" not in line:
        orderparams.append(line.split())
  else:
    for prm in argDict["ORDERPARAMS"].split(">"):
      orderparams.append(prm.split(","))
  if argDict["NODELIST"]=="read":
    nodelist = [line.split("\n")[0] for line in open("nodelist").readlines()]
  else:
    for node in argDict["NODELIST"].split(","):
      nodelist.append(node)
  for i in range(len(orderparams)-1):
    tmp = orderparams[i]
    dirname = currDir + "/" + fname + "-%03d-%03d"% ((int(float(tmp[0])*100)), (int(float(tmp[1])*100)))
    if os.path.isfile(os.path.join(dirname, argDict["TINKERXYZ"].split(".xyz")[0]+".bar")):
      print(" .bar file exists in %s !"%dirname)
      subprocess.run("rm -f %s"%os.path.join(dirname, argDict["TINKERXYZ"].split(".xyz")[0]+".bar"),shell=True)
      sub(nodelist[i], argDict["BARRUN"], dirname)
      print(" Deleted .bar file in %s and Resubmitted bar job on %s !"%(dirname, nodelist[i]))
    else:
      sub(nodelist[i], argDict["BARRUN"], dirname) 
      print(" Submitted bar job on %s !"%nodelist[i])
  return

def result():
  currDir = os.getcwd()
  argDict = readini()
  Free = []
  Err = []
  Dir = []
  orderparams = []
  fname = argDict["DIRNAME"]
  if argDict["ORDERPARAMS"] == "read":
    lines = open("orderparams").readlines()
    for line in lines:
      if "#" not in line:
        orderparams.append(line.split())
  else:
    for prm in argDict["ORDERPARAMS"].split(">"):
      orderparams.append(prm.split(","))
  for i in range(len(orderparams)-1):
    tmp = orderparams[i]
    dirname = currDir + "/" + fname + "-%03d-%03d"% ((int(float(tmp[0])*100)), (int(float(tmp[1])*100)))
    Dir.append(fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
    if os.path.isfile(os.path.join(dirname, "freeEnergy.log")):
      for line in open(os.path.join(dirname,"freeEnergy.log")).readlines():
        if "Free Energy via BAR Iteration" in line:
          Free.append(float(line.split()[-4]))
          Err.append(float(line.split()[-2]))
    else:
          print("Free Energy calculation not complete in %s"%fname+"-%03d-%03d"%((int(float(tmp[0])*100)), (int(float(tmp[1])*100))))
  for i,j,k in zip(Dir,Free,Err):
    print("%20s%10.2f%10.2f"%(i,j,k))
  totFree = np.array(Free).sum()
  totErr = ((np.array(Err)**2).sum())**0.5
  print("%20s%10.2f%10.2f"%("Total Free Energy", totFree, totErr))
  return
    
def main():
  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result}
  for action in sys.argv[1:]:
    if action in actions:
      actions[action]()

if __name__ == '__main__':
  main()
