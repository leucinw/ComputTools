
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import subprocess
import argparse
import numpy as np
import concurrent.futures 
from tqdm import tqdm
from datetime import datetime

'''check the availability of ONE node'''
def checkone(node):
  njobs = 0
  exelist = ["psi4","g09","g16","dynamic.x","cp2k.ssmp","mpirun_qchem","orca_mp2_mpi","gmx_mpi","gmx"]
  try:
    topstr = subprocess.check_output(f"ssh {node} 'top -n1 -b'", shell=True).decode("utf-8")
    line = ''.join(list(topstr))
    for exe in exelist:
      njobs += line.count(exe)
  except:
    pass
  return node,njobs

'''check the availability of ALL node'''
def checkNodes(nodes):
  now = datetime.now().strftime("%b %d %Y %H:%M:%S")
  print('\033[91m' + "[" + now + "] " + "checking availability of renlab clusters..." + '\033[0m')
  nodejobs = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(checkone, node) for node in nodes]
    for f in concurrent.futures.as_completed(results):
      nodejobs.append(f.result())
  return nodejobs 

'''prepare input filelist and nodelist that meets users requirements'''
def prepare(flistin, scratch=300, memory=30, maxmem=999):
  flist = []
  nlist = []
  nodes = np.loadtxt("/home/liuchw/shared/renlab.nodes", usecols=(0), unpack=True, dtype="str", skiprows=1) 
  memorys, cpus, scratches = np.loadtxt("/home/liuchw/shared/renlab.nodes", usecols=(1,2,3), unpack=True, dtype="int", skiprows=1)
  for f in flistin:
    if not os.path.isfile(os.path.splitext(f)[0] + ".log"):
      flist.append(f)
  for n, s, m in zip(nodes, scratches, memorys):
    if (s > scratch) and (m > memory) and (m < maxmem):
      nlist.append(n)
  return flist, nlist

'''submit ONE qm job to ONE node'''
def subOneJob(node,qmfile):
  f, ext = os.path.splitext(qmfile)
  cwd = os.getcwd()
  ext = ext.upper()
  if ext == ".COM":
    srcfile = "/home/liuchw/.bashrc.G09"
    exestr = "nohup g09 %s.com %s.log >log.sub 2>err.sub &"%(f, f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  elif ext == ".PSI4":
    srcfile = "/home/liuchw/.bashrc.psi4"
    exestr = "nohup psi4 -n 8 -i %s.psi4 -o %s.log >log.sub 2>err.sub &"%(f, f)
    cmdstr = 'ssh %s "source %s; cd %s; %s" &' % (node, srcfile, cwd, exestr)
  else:
    cmdstr = "echo 'file format %s not supported!'"%ext
  subprocess.run(cmdstr,shell=True)
  return


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'input',  nargs='+', help = "Input files",    required=True)  
  parser.add_argument('-n', dest = 'nodes',  nargs='+', help = "Submit jobs on these nodes ONLY. Default: []",default=[])  
  parser.add_argument('-x', dest = 'xnodes', nargs='+', help = "Submit jobs NOT on these nodes. Default: []", default=[])  
  parser.add_argument('-d', dest = 'disk',   type =int, help = "Disk of requested nodes. Default: 200 (GB)",  default=200)  
  parser.add_argument('-m', dest = 'memory', type =int, help = "Memory lower bound. Default: 30 (GB)", default=30)  
  parser.add_argument('-M', dest = 'maxmem', type =int, help = "Memory upper bound. Default: 999 (GB)", default=999)  
  args = vars(parser.parse_args())
  inps = args["input"]
  nodes = args["nodes"]
  xnodes = args["xnodes"]
  disk = args["disk"]
  memory = args["memory"] 
  maxmem = args["maxmem"] 

  # prepare flist, nlist
  if (nodes == []):
    qmfiles, nodes = prepare(inps, disk, memory, maxmem)
  # prepare flist only if nlist is specified
  else:
    qmfiles, _ = prepare(inps, disk, memory, maxmem)
  # remove the nodes user want to exclude
  if (xnodes != []):
    for ex in xnodes:
      if (ex in nodes):
        nodes.remove(ex)

  print('\033[92mSatisfied Nodes LIST:\033[0m')
  for node in nodes:
    print("%8s"%node)
  # submit jobs
  jobidx = 0
  if len(qmfiles) != 0:
    while(jobidx != len(qmfiles)):
      results = checkNodes(nodes)
      idlenodes = []
      for r in results:
        if (r[1]==0):
          idlenodes.append(r[0])
      if (idlenodes == []):
        time.sleep(30.0)
      else:
        for node in idlenodes:
          if (jobidx == len(qmfiles)):
            break
          subOneJob(node, qmfiles[jobidx])
          print("[%4s / %4s] "%(jobidx+1, len(qmfiles)) + '\033[92mSubmitted %s on %s!\033[0m'%(qmfiles[jobidx], node))
          jobidx += 1
        if (jobidx == len(qmfiles)):
          break
        else:
          time.sleep(30.0)
  else:
    print('\033[91m' + "Outputs already exist for your inputs!!" + '\033[0m')
  return

if len(sys.argv) == 1:
  print('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
else:
  main()
