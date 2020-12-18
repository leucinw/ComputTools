
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Assign new atom types to the specified key and xyz files
# Usage: python latomtyper.py xyz key atomIndex
 
import sys
import numpy as np

RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'


def main(xyz, key, idx):
  #dealing with xyz
  atomNumbers = np.loadtxt(xyz, usecols=(5), skiprows=1, dtype="int", unpack=True)
  idxAdd = idx - min(atomNumbers)
  lines = open(xyz).readlines()
  typesDict = {}
  with open(xyz + "_1", "w") as f:
    f.write(lines[0])
    for line in lines[1:]:
      s = line.split()
      formatted = "%6s%3s%12.6f%12.6f%12.6f%5s  "%(s[0], s[1], float(s[2]), float(s[3]), float(s[4]), str(int(s[5])+idxAdd))
      newline = formatted + "   ".join(s[6:]) + "\n"
      if int(s[5]) not in typesDict:
        typesDict[int(s[5])] = int(s[5]) + idxAdd 
      f.write(newline) 
  #dealing with key 
  lines = open(key).readlines()
  with open(key + "_1","w") as f:
    for line in lines:
      s = line.split()
      if "atom " in line:
        formatted = "%4s    %5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), s[3])
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "polarize " in line:
        if len(s) == 4: 
          formatted = "%8s%5s%8.4f%8.4f   "%(s[0], str(typesDict[int(s[1])]), float(s[2]), float(s[3]))
          grps = ' ' 
        else:
          grps = []
          formatted = "%8s%5s%8.4f%8.4f   "%(s[0], str(typesDict[int(s[1])]), float(s[2]), float(s[3]))
          for grp in s[4:]:
            if (int(float(grp)) in typesDict):
              grps.append(str(typesDict[int(grp)]))
            else:
              grps.append(str(grp))
        newline = formatted + "   ".join(grps) + "\n"
      elif "vdw " in line:
        # Hydrogens
        if len(s) == 5:
          formatted = "%3s     %5s%8.4f%8.4f%8.4f   "%(s[0], str(typesDict[int(s[1])]), float(s[2]), float(s[3]), float(s[4]))
        else:
          formatted = "%3s     %5s%8.4f%8.4f   "%(s[0], str(typesDict[int(s[1])]), float(s[2]), float(s[3]))
        newline = formatted + "\n"
      elif "bond " in line and ("#" not in line):
        formatted = "%4s    %5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]))
        newline = formatted + "   ".join(s[3:]) + "\n"
      elif "cflux-b " in line:
        formatted = "%7s    %5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]))
        newline = formatted + "   ".join(s[3:]) + "\n"
      elif "cp " in line:
        formatted = "%2s    %5s   "%(s[0], str(typesDict[int(s[1])]))
        newline = formatted + "   ".join(s[2:]) + "\n"
      elif "ct " in line:
        formatted = "%2s    %5s   "%(s[0], str(typesDict[int(s[1])]))
        newline = formatted + "   ".join(s[2:]) + "\n"
      elif "opbend " in line:
        formatted = "%4s    %5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]))
        newline = formatted + "   ".join(s[3:]) + "\n"
      elif "strbnd " in line:
        formatted = "%5s    %5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), str(typesDict[int(s[3])]))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "angle " in line and ("#" not in line):
        formatted = "%5s    %5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), str(typesDict[int(s[3])]))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "cflux-a " in line:
        formatted = "%7s    %5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), str(typesDict[int(s[3])]))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "anglep " in line:
        formatted = "%6s    %5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), str(typesDict[int(s[3])]))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "torsion " in line:
        formatted = "%7s    %5s%5s%5s%5s   "%(s[0], str(typesDict[int(s[1])]), str(typesDict[int(s[2])]), str(typesDict[int(s[3])]), str(typesDict[int(s[4])]))
        newline = formatted + "   ".join(s[5:]) + "\n"
      elif "multipole " in line:
        mpoleframes = s[1:-1]
        for i in range(len(mpoleframes)):
          if int(mpoleframes[i]) == 0:
            mpoleframes[i] = 0 
          elif int(mpoleframes[i]) < 0:
            mpoleframes[i] = "-"+str(typesDict[abs(int(mpoleframes[i]))]) 
          else:
            mpoleframes[i] = str(typesDict[int(mpoleframes[i])]) 
        if len(mpoleframes) == 4:# Z-Bisector or 3-Fold
          formatted = "%9s    %5s%5s%5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1], mpoleframes[2], mpoleframes[3])
        elif len(mpoleframes) == 3:# Z-then-X or Bisector
          formatted = "%9s    %5s%5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1], mpoleframes[2])
        elif len(mpoleframes) == 2:# Z-only
          formatted = "%9s    %5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1])
        elif len(mpoleframes) == 2:# Charge only 
          formatted = "%9s    %5s   "%(s[0], mpoleframes[0])
        else:
          print("multipole" + mpoleframes)
        newline = formatted + "         " + s[-1] + "\n"
      else:
        newline = line
      f.write(newline) 
  return

# Execute
if len(sys.argv) != 4:
  sys.exit("\n" + RED + "   Usage: python latomtyper.py xyz key index" + ENDC + "\n")
else:
  xyz = sys.argv[1]
  key = sys.argv[2]
  idx = int(sys.argv[3])
  main(xyz, key, idx)
