
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import numpy as np

def readINPUTS():
  types = np.loadtxt(xyz, usecols=(5), dtype='str', unpack=True, skiprows=1)
  settypes = list(set(types))
  oldtypes = []
  oldclasses = []
  parameters = open(prm).readlines()
  for line in parameters:
    d = line.split()
    for ot in settypes:
      if ("atom " in line) and (d[1] == ot):
        oldtypes.append(d[1])
        oldclasses.append(d[2])
  newtypes = [str(i) for i in range(idx, idx+len(oldtypes))]
  newclasses = [str(i) for i in range(idx, idx+len(oldclasses))]
  type_maps = dict(zip(oldtypes, newtypes))
  class_maps = dict(zip(oldclasses, newclasses))
  
  lines = open(xyz).readlines()
  with open(newxyz, "w") as f:
    f.write(lines[0])
    for line in lines[1:]:
      d = line.split()
      newline = "%6s%3s%12.6f%12.6f%12.6f%6s  %s\n"%(d[0], d[1], float(d[2]), float(d[3]), float(d[4]), type_maps[d[5]], ' '.join(d[6:]))
      f.write(newline)
  return type_maps,class_maps,parameters

def ClassBasedTerms(pline, term):
  d = pline.split()
  with open(newprm, 'a') as f:
    if term == 'atom':
      if set([d[1]]).issubset(set(type_maps.keys())):
        newtypesclasses = [type_maps[d[1]], class_maps[d[2]]]
        newline = '  '.join([d[0], newtypesclasses[0], newtypesclasses[1], '  '.join(d[3:])])
        f.write(newline + "\n")
    elif term == 'bond':
      if set([d[1], d[2]]).issubset(set(class_maps.keys())):
        newclasses = [class_maps[d[1]], class_maps[d[2]]]
        newline = '  '.join([d[0], newclasses[0], newclasses[1], d[3], d[4]])
        f.write(newline + "\n")
    elif (term == "angle") or (term == "anglep"):
      if set([d[1], d[2], d[3]]).issubset(set(class_maps.keys())):
        newclasses = [class_maps[d[1]], class_maps[d[2]], class_maps[d[3]]]
        newline = '  '.join([d[0], newclasses[0], newclasses[1], newclasses[2], d[4], d[5]])
        f.write(newline + "\n")
    elif term == "strbnd":
      if set([d[1], d[2], d[3]]).issubset(set(class_maps.keys())):
        newclasses = [class_maps[d[1]], class_maps[d[2]], class_maps[d[3]]]
        newline = '  '.join([d[0], newclasses[0], newclasses[1], newclasses[2], d[4], d[5]])
        f.write(newline + "\n")
    elif term == 'opbend':
      if set([d[1], d[2]]).issubset(set(class_maps.keys())):
        newclasses = [class_maps[d[1]], class_maps[d[2]]]
        newline = '  '.join([d[0], newclasses[0], newclasses[1], d[3], d[4], d[5]])
        f.write(newline + "\n")
    elif term == "torsion":
      if set([d[1], d[2], d[3], d[4]]).issubset(set(class_maps.keys())):
        newclasses = [class_maps[d[1]], class_maps[d[2]], class_maps[d[3]], class_maps[d[4]]]
        newline = '  '.join([d[0], newclasses[0], newclasses[1], newclasses[2], newclasses[3]] + d[5:])
        f.write(newline + "\n")
    elif term == 'vdw':
      if d[1] in class_maps.keys():
        newclass = class_maps[d[1]]
        newline = '  '.join([d[0], newclass, '  '.join(d[2:])])
        f.write(newline + "\n")
    else:
      sys.exit(f"Error: {term} not supported by ClassBasedTerms()") 
  return

def TypeBasedTerms(pline, term):
  d = pline.split()
  groups = []
  with open(newprm, "a") as f: 
    if term == 'polarize':
      if d[1] in type_maps.keys():
        newtypes = type_maps[d[1]]
        if fun == "AMOEBA":
          grpidx = 3
        if fun == "AMOEBAPLUS":
          grpidx = 4
        for i in range(len(d)-1,grpidx,-1):
          if d[i] in type_maps.keys():
            groups.append(type_maps[d[i]])
          else:
            groups.append('0')
        groups.reverse()
        if groups != []:
          if fun == "AMOEBA":
            newline = '  '.join([d[0], newtypes, '  '.join(d[2:4])] + groups)
          if fun == "AMOEBAPLUS":
            newline = '  '.join([d[0], newtypes, '  '.join(d[2:5])] + groups)
        else:
          if fun == "AMOEBA":
            newline = '  '.join([d[0], newtypes, '  '.join(d[2:4])])
          if fun == "AMOEBAPLUS":
            newline = '  '.join([d[0], newtypes, '  '.join(d[2:5])])
        f.write(newline + "\n")
    
    elif term == "multipole":
      index = parameters.index(pline)
      newtypes = []
      if d[1] in type_maps.keys():
        newtypes.append(type_maps[d[1]])
        for i in range(2, len(d)-1):
          if int(d[i]) > 0:
            if d[i] in type_maps.keys():
              newtypes.append(type_maps[d[i]])
            else:
              newtypes.append('0')
          elif int(d[i]) < 0:
            if str(abs(int(d[i]))) in type_maps.keys():
              newtypes.append("-" + type_maps[str(abs(int(d[i])))])
            else:
              newtypes.append('-0')
          else:
            newtypes.append('0')
        newline = '         '.join([d[0], '  '.join(newtypes), d[-1]])
        f.write(newline + "\n")
        for i in range(index+1, index+5):
          f.write(parameters[i])
    else:
      sys.exit(f"Error: {term} not supported by TypeBasedTerms()")
  
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz', dest = 'xyz', required=True, help="tinker xyz")
  parser.add_argument('-prm', dest = 'prm', required=True, help="tinker key")
  parser.add_argument('-idx', dest = 'idx', required=True, help="atom index", type=int)
  parser.add_argument('-fun', dest = 'fun', default='AMOEBAPLUS', help="force field", type=str.upper, choices=['AMOEBA', 'AMOEBAPLUS'])
  args = vars(parser.parse_args())
  global xyz,prm,idx,fun 
  xyz = args["xyz"]
  prm = args["prm"]
  idx = args["idx"]
  fun = args["fun"]

  global type_maps
  global class_maps
  global parameters
  global newxyz,newprm

  newxyz = xyz + "_1"
  newprm = prm + "_1"
  if os.path.isfile(newxyz):
    os.remove(newxyz)
  if os.path.isfile(newprm):
    os.remove(newprm)
  
  databasedir="/home/liuchw/shared/forcefield/" 
  if fun == "AMOEBA":
    head = open(os.path.join(databasedir, "amoebabio18_header.prm")).readlines()
  else:
    head = open(os.path.join(databasedir, "amoebaplus21_header.prm")).readlines()
  with open(newprm, "w") as f:
    for h in head:
      f.write(h)

  type_maps, class_maps, parameters = readINPUTS()
  keywords = ["atom ", "bond ", "angle ", "anglep ", "strbnd ", "opbend ", "torsion ", "vdw ", "polarize ", "multipole "]
  valwords = [k.split()[0] for k in keywords]
  keyvals = dict(zip(keywords,valwords))
  classBased_terms = dict(list(keyvals.items())[:-2]) 
  typeBased_terms = dict(list(keyvals.items())[-2:])

  maxtype = max(type_maps.values())
  print(idx, maxtype)
  for p in parameters:
    for t in classBased_terms: 
      if t in p:
        ClassBasedTerms(p, classBased_terms[t])
    for t in typeBased_terms: 
      if t in p:
        TypeBasedTerms(p, typeBased_terms[t])
  return

if __name__ == "__main__":
  main()