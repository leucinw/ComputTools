
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import os,time,sys,subprocess
from colorama import Fore
from multiprocessing import Pool

usage = ''' Currently supported syntax:
          1. python lsub.py water1.com   [water2.com]   [..]
          2. python lsub.py water1.qchem [water2.qchem] [..]
          3. python lsub.py water1.psi4  [water2.psi4]  [..]
          4. python lsub.py *.com
          5. python lsub.py *.qchem
          6. python lsub.py *.psi4
          7. python lsub.py dynamic_omm.cuda
  '''

# global list 
CPU_ExeList = ["psi4","g09","g16","dynamic.x","cp2k.ssmp","mpirun_qchem", "orca_mp2_mpi"]
GPU_ExeList = ["dynamic_omm.x","bar_omm.x", "pmemd"]

def checkOneNode(eachNode, GPU_Job):
  nCPUjobs = 0
  nGPUjobs = 0
  cmdstr = "ssh %s 'top -n1 -b'"%eachNode
  try:
    topstr = subprocess.check_output(cmdstr, shell=True).decode("utf-8")
    linestr = ''.join(list(topstr))
    for exe in GPU_ExeList:
      nGPUjobs += linestr.count(exe)
    for exe in CPU_ExeList:
      nCPUjobs += linestr.count(exe)
  except:
    pass
  if GPU_Job:
    return eachNode, nGPUjobs
  else:
    return eachNode, nCPUjobs

def checkAllNodes(nodelistfile, GPU_Job):
  print(Fore.RED + "Checking nodes ...")
  Nodes = []; nJobs = []; idleNodes = []
  for line in open(nodelistfile).readlines():
    if "#" not in line[:1]:
      d = line.split()
      Nodes.append([d[0], GPU_Job])
      nJobs.append(int(d[1]))
  p = Pool(len(Nodes))
  results = dict(p.starmap(checkOneNode,Nodes))
  p.close()
  p.join()
  for [node,_], nJob in zip(Nodes, nJobs):
    num = nJob - results[node] 
    idleNodes += [node]*num
  return idleNodes

def checkMem(node, nodelistfile):
  for line in open(nodelistfile).readlines():
    if ("#" not in line) and (line !="\n"):
      if node in line:
        mem = int(line.split()[2])
  return mem


def checkJobs(filelist):
  fileToSubmit = []
  fileformats = [".qchem", ".psi4", ".com"]
  for files in filelist:
    _, ext = os.path.splitext(files)
    if ext in fileformats: 
      fname, end = os.path.splitext(files)
      if not os.path.isfile(fname + ".log"):
        fileToSubmit.append(files)
  return fileToSubmit 

def subOneJob(node, inputFile):
  f, ext = os.path.splitext(inputFile)
  cwd = os.getcwd()
  # qchem 
  if ext == ".qchem":
    srcfile = "/home/liuchw/.bashrc.qchem"
    exestr = "nohup qchem -nt 8 %s.qchem %s.log >log.sub 2>err.sub &"%(f,f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &'%(node,srcfile,cwd,exestr)
  # psi4 
  elif ext == ".psi4":
    srcfile = "/home/liuchw/.bashrc.poltype"
    jobtype = "psi4"
    memmax = checkMem(node, "/home/liuchw/bin/Psi4Node")
    lines = open(inputFile).readlines()
    for line in lines:
      if "MEM" in line.upper():
        mem = int(line.split()[1])
    if mem < memmax:
      exestr = "nohup psi4 -n 8 -i %s.psi4 -o %s.log >log.sub 2>err.sub &"%(f, f)
      cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
    else:
      cmdstr = "echo 'memory not big enough on %s !' "%node
  # gaussian
  elif ext == ".com":
    srcfile = "/home/liuchw/.bashrc.G09"
    exestr = "nohup g09 %s.com %s.log >log.sub 2>err.sub &"%(f, f)
    #srcfile = "/home/liuchw/.bashrc.G16"
    #exestr = "nohup g16 %s.com %s.log >log.sub 2>err.sub &"%(f, f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  # cuda 
  elif ext == ".cuda":
    exestr = "nohup sh %s >log.sub 2>err.sub &"%inputFile
    cmdstr = 'ssh %s "cd %s; %s" &' % (node, cwd, exestr)
  else:
    cmdstr = "echo 'file format not supported!'"
  subprocess.run(cmdstr,shell=True)
  return

def subMultipleJobs(filelist, nodelistfile, GPU_Job):
  currentJob = 0
  checkTime = 0.5
  if not GPU_Job: 
    files = checkJobs(filelist)
  else:
    files = filelist
  while(currentJob != len(files)):
    idleNodes = checkAllNodes(nodelistfile, GPU_Job)
    if (idleNodes == []):
      #print(Fore.RED + "No idle nodes ! Wait for %s minute(s) !"%checkTime)
      time.sleep(60.0*checkTime)
    else:
      for eachNode in idleNodes:
        if (currentJob == len(files)):break
        subOneJob(eachNode,files[currentJob])
        print(Fore.GREEN + "Submitted %s on %s!"%(files[currentJob], eachNode))
        currentJob += 1
      if (currentJob == len(files)):
        break
      else:
        time.sleep(60.0*checkTime)
  return 

def main():
  filelist = []
  jobtype = ''
  if len(sys.argv) > 1:
    filelist = sys.argv[1:]
    _, ext = os.path.splitext(filelist[0])
    jobtype = ext[1:].upper()
  else:
    print(Fore.RED + usage)

  if jobtype == "CUDA":
    GPU_Job = True
  else:
    GPU_Job = False

  supportedTypes = ["COM", "PSI4",  "QCHEM", "CUDA"]
  nodeDict = {"COM"  : "/home/liuchw/bin/GaussianNode",\
              "PSI4" : "/home/liuchw/bin/Psi4Node", \
              "QCHEM": "/home/liuchw/bin/QChemNode", \
              "CUDA":  "/home/liuchw/bin/CudaNode"}
  if (jobtype not in supportedTypes):
    print(Fore.RED + "%s files not supported!" %jobtype)
    sys.exit(1)
  if (filelist == []):
    print(Fore.RED + "you may have broken .log files?")
  if (filelist != []) and (jobtype in supportedTypes):
    pairs = []
    for f in filelist:
      size = os.path.getsize(f)
      pairs.append((size,f))
    pairs.sort(key=lambda s: s[0])
    sortedfilelist=[]
    for pair in pairs:
      sortedfilelist.append(pair[1])
    subMultipleJobs(sortedfilelist, nodeDict[jobtype], GPU_Job)
  return 

if __name__ == '__main__':
  main()
