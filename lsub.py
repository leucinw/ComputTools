
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import os,time,sys,subprocess
from colorama import Fore

usage = ''' Currently supported syntax:
  1. python lsub.py water1.com   [water2.com]   [..]
  2. python lsub.py water1.qchem [water2.qchem] [..]
  3. python lsub.py water1.psi4  [water2.psi4]  [..]
  4. python lsub.py *.com
  5. python lsub.py *.qchem
  6. python lsub.py *.psi4
  '''

def checkNode(nodelistfile, checkExeList):
  Nodes = []; nJobs = []; idleNodes = []
  for line in open(nodelistfile).readlines():
    if "#" not in line[:1]:
      d = line.split()
      Nodes.append(d[0])
      nJobs.append(int(d[1]))
  for node,nJob in zip(Nodes, nJobs):
    njobs = 0
    subprocess.run("rm -rf top.tmp", shell=True)
    cmdstr = "ssh %s 'top -n1 -b' >top.tmp"%node
    subprocess.run(cmdstr,shell=True)
    linestr = ''.join(open("top.tmp").readlines())
    for exe in checkExeList:
      njobs += linestr.count(exe)
    num = nJob - njobs
    idleNodes += [node]*num
  return idleNodes

def checkJobs(filelist):
  fileToSubmit = []
  fileformats = [".qchem", ".psi4", ".com"]
  _, ext = os.path.splitext(open(filelist).readlines()[0][:-1])
  if ext in fileformats: 
    for f in open(filelist).readlines():
      fname, end = os.path.splitext(f[:-1])
      if not os.path.isfile(fname + ".log"):
        fileToSubmit.append(f[:-1])
  return fileToSubmit 

def subOneJob(node, input_file):
  f, ext = os.path.splitext(input_file)
  cwd = os.getcwd()
  # qchem 
  if ext == ".qchem":
    srcfile = "/home/liuchw/.bashrc.qchem"
    exestr = "nohup qchem -nt 8 %s.qchem %s.log >log.sub 2>err.sub &"%(f,f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &'%(node,srcfile,cwd,exestr)
  # psi4 
  if ext == ".psi4":
    srcfile = "/home/liuchw/.bashrc.poltype"
    jobtype = "psi4"
    exestr = "nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.sub &"%(f, f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  # gaussian
  if ext == ".com":
    srcfile = "/home/liuchw/.bashrc.G16"
    exestr = "nohup g16 %s.com %s.log 2>err.sub &"%(f, f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  subprocess.run(cmdstr,shell=True)
  return

def subMultipleJobs(filelist, nodelistfile, checkExeList):
  currentJob = 0
  checkTime = 1.0 
  files = checkJobs(filelist)
  while(currentJob != len(files)):
    idleNodes = checkNode(nodelistfile, checkExeList)
    if (idleNodes == []):
      print(Fore.RED + "No idle nodes ! Wait for %s minute(s) !"%checkTime)
      time.sleep(60.0*checkTime)
    else:
      for eachNode in idleNodes:
        if (currentJob == len(files)):break
        subOneJob(eachNode,files[currentJob])
        print(Fore.GREEN + "Submitted %s on %s!"%(files[currentJob], eachNode))
        currentJob += 1
      time.sleep(60.0*checkTime)
  return 

def main():
  filelist = []
  jobtype = ''
  if len(sys.argv) > 1:
    if len(sys.argv) == 2:
      fname, ext = os.path.splitext(sys.argv[1])
      jobtype = ext[1:].upper()
      if fname == "*":
        files = os.listdir(os.getcwd())
        for eachfile in files:
          if os.path.splitext(eachfile)[-1] == ext:
            filelist.append(eachfile)
      else:
        filelist.append(sys.argv[2])
    else:
      filelist = sys.argv[2:]
  else:
    print(Fore.RED + usage)
 

  checkExeList = ["psi4","g09","g16","dynamic.x","cp2k.ssmp","mpirun_qchem"]
  nodeDict = {"COM"  : "/home/liuchw/bin/GaussianNode",\
              "PSI4" : "/home/liuchw/bin/Psi4Node", \
              "QCHEM": "/home/liuchw/bin/QChemNode"}

  if (filelist  != []):
    subMultipleJobs(filelist, nodeDict[jobtype], checkExeList)
  return 

if __name__ == '__main__':
  main()
