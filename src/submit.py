#!/usr/bin/env python

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import subprocess
import concurrent.futures 

RED = '\033[91m'
ENDC = '\033[0m'

''' submit one or more sh files to one node '''

def subone(script):
  subcmd = f"cd {path}; sh {script}"
  os.system(subcmd)
  return True

def submore(scripts):
  jobstats = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(subone, script) for script in scripts]
    for f in concurrent.futures.as_completed(results):
      jobstats.append(f.result())
  return jobstats

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', dest = 'path', help = "path to submit", default = ' ')  
  parser.add_argument('-x', dest = 'scripts', nargs='+', help = "executable script files", default = [], required = True)  
  args = vars(parser.parse_args())
  global path,scripts
  path = args["path"]
  scripts = args["scripts"]
  
  r = subprocess.check_output("hostname", shell=True).decode("utf-8")
  if ('bme-nova.bme.utexas.edu' in r):
    sys.exit(RED + "Please run this on a computing node" + ENDC)

  if path == ' ':
    path = os.getcwd()
  else:
    path = os.path.join(os.getcwd(), path)

  if len(scripts) == 1:
    subone(scripts[0])
  else:
    submore(scripts)
  return

if __name__ == "__main__":
  main()
