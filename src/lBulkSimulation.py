# write bulk(and gas phase) simulations to obtain
# 1. Density
# 2. Hvap

# files to prep:
# 1. liquid.xyz
# 2. gas.xyz
# 3. gas.key
# 4. liquid.key
# 5. amoeba09.prm

import sys
import subprocess
import os

# check if files exist
def checkfiles():
  MODE = 1
  checklist = ["liquid.key", "liquid.xyz", "gas.xyz", "gas.key"]
  files = os.listdir(os.getcwd())
  if checklist[0] in files:
    lines = open(checklist[0]).readlines()
    for line in lines:
      if "parameters" in line.lower():
        prmfile = line.split()[1][:-1]
        break
    checklist.append(prmfile)
  # check MODE 1
  for f in checklist:
    if f not in files:
      MODE = 0
  # check MODE 2
  if MODE == 1:
    checklist = ["liquid.arc", "liquid.log", "gas.arc", "gas.log"]
    files = os.listdir(os.getcwd())
    for f in checklist:
      if f not in files:
        MODE = 0
    MODE = 2
  

  return MODE

# run gas phase dynamic
def gas_dynamic():

# run gas phase analyze
def gas_analyze():

# run liquid dynamic
def liquid_dynamic():

# run liquid analyze
def liquid_analyze():

# calc liquid density from .arc
def liquid_density():

# calc liquid hvap 
def hvap():

# main 
def main():
  
# execute
if __main__ == "__main__":
  main()
