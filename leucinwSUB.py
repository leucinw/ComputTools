#!/usr/bin/env python3

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os,sys
import subprocess

def sub_one_job(infile, node=None,userpath=None):
    data = os.path.splitext(infile)
    fname = data[0]
    ext = data[1]
    if userpath == None: 
      cwd = os.getcwd()
    else:
      cwd = userpath
    # psi4
    if (ext == ".psi4") or (ext == ".psi"):
      src = "/home/liuchw/.bashrc.psi4conda"
      if node==None:
        cmdstr = "source %s && cd %s && nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.log &" %(src,cwd,fname,fname)
      else:
        cmdstr = 'ssh %s "source %s; cd %s; nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.log &" &'%(node,src,cwd,fname,fname)
    # qchem 
    if ext == ".qchem":
      src = "/home/liuchw/.bashrc.qchem"
      if node==None:
        cmdstr = "source %s; cd %s; nohup qchem -np 8 %s.qchem %s.log 2>err.log &" %(src,cwd,fname,fname)
      else:
        cmdstr = 'ssh %s "source %s; cd %s; nohup qchem -np 8 %s.qchem %s.log 2>err.log &" &'%(node,src,cwd,fname,fname)
    # gaussian
    if (ext == ".com") or (ext == ".gjf"):
      src = "/home/liuchw/.bashrc.G09"
      if node == None:
        cmdstr = "source %s; cd %s; nohup g09 %s.com %s.out 2>err.log &" % (src, cwd, fname, fname)
      else:
        cmdstr = 'ssh %s "source %s; cd %s; nohup g09 %s.com %s.out 2>err.log &" &' % (node, src, cwd, fname, fname)
    # poltype 
    if ext == ".poltype":
      src = "/home/liuchw/.bashrc.poltype"
      if node == None:
        cmdstr = "source %s; cd %s; nohup sh %s.poltype 2>err.log &"%(node, src, cwd, fname)
      else:
        cmdstr = 'ssh %s "source %s; cd %s; nohup sh %s.poltype 2>err.log &" &' % (node, src, cwd, fname)

    # normal executable .sh
    if ext == ".sh":
      if node == None:
        cmdstr = "cd %s; nohup sh %s.sh 2>err.log &" % (cwd, fname)
      else:
        cmdstr = 'ssh %s "cd %s; nohup sh %s.sh 2>err.log &" &' % (node, cwd, fname)
    subprocess.run(cmdstr,shell=True)
    return

def main():
  usage = '''        ===============================================================
        qsub.py is a python script to submit jobs on Ren-Lab clusters.
        Files supported: psi4, gaussian, qchem, poltype.
        $ leucinwSUB.py file.ext [node] [path]
        where file.ext is the input file, node and path are optional.
        ==============================================================='''

  if len(sys.argv) == 2:
    sub_one_job(sys.argv[1])
    print("Submitted %s on the current node !"%(sys.argv[1]))
  elif len(sys.argv) == 3:
    sub_one_job(sys.argv[1], sys.argv[2])
    print("Submitted %s on %s !"%(sys.argv[1], sys.argv[2]))
  elif len(sys.argv) == 4:
    sub_one_job(sys.argv[1], sys.argv[2], sys.argv[3])
    print("Submitted %s on %s !"%(sys.argv[1], sys.argv[2]))
  else:
    print(usage)

if __name__ == "__main__":
  main()
