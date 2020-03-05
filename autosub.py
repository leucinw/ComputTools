import os,time,sys,subprocess

def checkNode():
  Nodes = []; nJobs = []; maxMem = []; idleNodes = []
  checkExeList = ["psi4", "g09", "g16", "dynamic.x", "cp2k.ssmp"]
  for line in open("/home/liuchw/bigMemoryNode").readlines():
    if "#" not in line[:1]:
      d = line.split()
      Nodes.append(d[0])
      nJobs.append(int(d[1]))
      maxMem.append(int(d[2]))
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
  return idleNodes,maxMem

def subOneJob(node, memory, input_file):
  fname, ext = os.path.splitext(input_file)
  cwd = os.getcwd()
  # qchem 
  if ext == ".qchem":
    srcfile = "/home/liuchw/.bashrc.qchem"
    exestr = "nohup qchem -np 8 %s.qchem %s.log 2>err.sub %"%(fname,fname)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &'%(node,src,cwd,exestr)
  # psi4 
  if ext == ".psi4":
    srcfile = "/home/liuchw/.bashrc.poltype"
    jobtype = "psi4"
    exestr = "nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.log &"%(fname, fname)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  # gaussian
  if ext == ".com":
    tempfile = input_file + "_1"
    f = open(tempfile, 'w')
    for line in open(input_file).readlines():
      if "MEM" in line.upper():
        if int(line.split("=")[-1][:-3]) > memory:
          allowedMem = int(memory*0.5)
          line = "%Mem=" + "%sGB\n"%(allowedMem)
      f.write(line)
    f.close()
    os.rename(tempfile, input_file)
    srcfile = "/home/liuchw/.bashrc.G16"
    exestr = "nohup g16 %s.com %s.out 2>err.log &"%(fname, fname)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  subprocess.run(cmdstr,shell=True)
  return

def subMultipleJobs(filelist, checkTime):
  currentJob = 0
  files = [line[:-1] for line in open(filelist).readlines()]
  while(True):
    idleNodes = []
    if currentJob == (len(files)):
      break
    idleNodes, memoryList = checkNode()
    if idleNodes == []:
      if checkTime > 1:
        print("I will sleep for %s minutes !"%checkTime)
      else:
        print("I will sleep for %s minute !"%checkTime)
      time.sleep(60.0*checkTime)
    else: 
      for eachNode, memory in zip(idleNodes, memoryList):
        if currentJob == (len(files)):
          break
        subOneJob(eachNode, memory, files[currentJob])
        print("I submitted %s on %s!"%(files[currentJob], eachNode))
        currentJob += 1
  return 

subMultipleJobs(sys.argv[1], float(sys.argv[2]))
