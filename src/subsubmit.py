
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import argparse

''' submit one or more sh files to one or more node(s).
    submit.py is required to run this program
'''

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-x', dest = 'scripts', nargs='+', help = "Excutable script. One or Multiple. ", required=True)  
  parser.add_argument('-n', dest = 'nodes', nargs='+', help = "Node list. One or Multipole. ",default=[])  
  args = vars(parser.parse_args())
  global nodes,path,scripts
  scripts = args["scripts"]
  nodes = args["nodes"]
  path = os.getcwd()
  
  nnode = len(nodes)
  nscript = len(scripts)
  nfile = int(nscript/nnode)
  os.system("rm -rf node*.sh")
  ss = [scripts[i*nscript // nnode: (i+1)*nscript // nnode] for i in range(nnode)]
  for i in range(nnode):
    nodestr = ' '.join(ss[i])
    with open(f"{nodes[i]}.sh", "w") as f:
      f.write(f'ssh {nodes[i]} "cd {path}; submit -x {nodestr}"\n')
  os.system("parallel sh ::: node*.sh")  
  return

if __name__ == "__main__":
  main()
