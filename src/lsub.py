
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os,time,sys,subprocess
import argparse
import numpy as np
import concurrent.futures 

# nodes to check
def checkone(node):
  njobs = 0
  exelist = ["psi4","g09","g16","dynamic.x","dynamic", "cp2k.ssmp","mpirun_qchem","orca_mp2_mpi","gmx_mpi","gmx"]
  try:
    topstr = subprocess.check_output(f"ssh {node} 'top -n1 -b'", shell=True).decode("utf-8")
    line = ''.join(list(topstr))
    for exe in exelist:
      njobs += linestr.count(exe)
  except:
    pass
  return node,njobs

def checkNodes(nodes):
  nodejobs = []
  print('\033[91m' + "Checking availability of RENLAB clusters..." + '\033[0m')
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(checkone, node) for node in nodes]
    for f in concurrent.futures.as_completed(results):
      nodejobs.append(f.result())
  return nodejobs 

# submit job one by one 
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

# prepare file&nodelist
def prepare(flistin, scratch=300, memory=30):
  flist = []; nlist = []
  nodes = np.loadtxt("/home/liuchw/shared/renlab.nodes", usecols=(0), unpack=True, dtype="str", skiprows=1) 
  memorys, cpus, scratches = np.loadtxt("/home/liuchw/shared/renlab.nodes", usecols=(1,2,3), unpack=True, dtype="int", skiprows=1)
  #check if .log file exists!
  for f in flistin:
    if not os.path.isfile(os.path.splitext(f)[0] + ".log"):
      flist.append(f)
    else:
      print('\033[93m' + f"log file exists for {f}!" + '\033[0m') 
  #select nodes according to scratch/memory
  for n, s, m in zip(nodes, scratches, memorys):
    if (s > scratch) and (m > memory):
      nlist.append(n)
  return flist, nlist

# main
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'input',  nargs='+', help = "Input files [x.com, y.psi4, z.qchem]",    required=True)  
  parser.add_argument('-n', dest = 'nodes',  nargs='+', help = "Submit jobs on these nodes ONLY. Default: []",default=[])  
  parser.add_argument('-x', dest = 'nodes_', nargs='+', help = "Submit jobs NOT on these nodes. Default: []", default=[])  
  parser.add_argument('-d', dest = 'disk',   type =int, help = "Disk of requested nodes. Default: 200 (GB)",  default=200)  
  parser.add_argument('-m', dest = 'memory', type =int, help = "Memory of requested nodes. Default: 30 (GB)", default=30)  
  args = vars(parser.parse_args())
  inps = args["input"]
  nodes = args["nodes"]
  nodes_ = args["nodes_"]
  disk = args["disk"]
  memory = args["memory"] 

  # use nodes
  if (nodes == []):
    qmfiles, nodes = prepare(inps, disk, memory)
  else:
    qmfiles, _ = prepare(inps, disk, memory)

  # excluding nodes
  if (nodes_ != []):
    for ex in nodes_:
      if (ex in nodes):
        nodes.remove(ex)

  # detect and submit jobs 
  jobidx = 0
  while(jobidx != len(qmfiles)):
    results = checkNodes(nodes)
    idlenodes = []
    for r in results:
      if (r[1]==0):
        idlenodes.append(r[0])
    if (idlenodes == []):time.sleep(30.0)
    else:
      for node in idlenodes:
        if (jobidx == len(qmfiles)):break
        subOneJob(node,qmfiles[jobidx])
        print('\033[92m' + "Submitted %s on %s!"%(qmfiles[jobidx],node) + '\033[0m')
        jobidx += 1
      if (jobidx == len(qmfiles)):break
      else:time.sleep(30.0)
  return

# execute
if __name__ == '__main__':
  if len(sys.argv) == 1:
    sys.exit('\033[93m' + " please use '-h' option to see usage" + '\033[0m')
  main()
