import os
import os.path
import sys


def subJobOnNode(node, input_file, userpath=None):
    data = os.path.splitext(input_file)
    fname = data[0]
    ext = data[1]
    if userpath == None: 
        cwd = os.getcwd()
    if userpath != None:
        cwd = userpath 
    if ext == ".psi4": # psi4 
        cmdstr = 'ssh %s "source /home/liuchw/.bashrc.psi4conda; cd %s; nohup psi4 -n 8 -i %s.psi4 -o %s.log 2>err.log &" &' % (
            node, cwd, fname, fname)
    if ext == ".com": # gaussian 
        cmdstr = 'ssh %s "source /home/liuchw/.bashrc.G16; cd %s; nohup g16 %s.com %s.out 2>err.log &" &' % (
            node, cwd, fname, fname)
    if ext == ".poltype": # poltype
        cmdstr = 'ssh %s "source ~/.bashrc.poltype; cd %s; nohup sh %s.sh 2>err.log &" &' % (
            node, cwd, fname)
    if ext == ".sh": # regular  
        cmdstr = 'ssh %s "cd %s; nohup sh %s.sh 2>err.log &" &' % (node, cwd, fname)
    os.system(cmdstr)
    return
