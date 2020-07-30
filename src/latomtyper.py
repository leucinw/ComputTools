
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

# Assign new atom types to the specified key and xyz files
# Usage: python latomtyper.py xyz key atomIndex
 
import sys
import numpy as np

xyz = sys.argv[1]
key = sys.argv[2]
idx = int(sys.argv[3])

def modXYZ(xyz):
  atomNumbers = np.loadtxt(xyz, usecols=(5), skiprows=1, dtype="int", unpack=True)
  idxAdd = idx - min(atomNumbers)
  lines = open(xyz).readlines()
  with open("ttt.xyz", "w") as f:
    f.write(lines[0])
    for line in lines[1:]:
      s = line.split()
      formatted = "%6s%3s%12.6f%12.6f%12.6f%5s  "%(s[0], s[1], float(s[2]), float(s[3]), float(s[4]), str(int(s[5])+idxAdd)) 
      newline = formatted + "   ".join(s[6:]) + "\n"
      f.write(newline) 
  return idxAdd

def modKEY(key):
  idxAdd = modXYZ(xyz)
  lines = open(key).readlines()
  with open("ttt.key","w") as f:
    for line in lines:
      s = line.split()
      if "atom " in line:
        formatted = "%4s    %5s%5s%5s   "%(s[0], str(int(s[1])+idxAdd), str(int(s[2])+idxAdd), s[3])
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "polarize " in line:
        formatted = "%8s%5s%8.4f%8.4f   "%(s[0], str(int(s[1])+idxAdd), float(s[2]), float(s[3]))
        # Mono-atomic ion
        if len(s) == 4:
          newline = formatted + "\n"
        else:
          grps = [str(int(grp)+idxAdd) for grp in s[4:]]
          newline = formatted + "   ".join(grps) + "\n"
      elif "vdw " in line:
        # Hydrogens
        if len(s) == 5:
          formatted = "%3s     %5s%8.4f%8.4f%8.4f   "%(s[0], str(int(s[1])+idxAdd), float(s[2]), float(s[3]), float(s[4]))
        else:
          formatted = "%3s     %5s%8.4f%8.4f   "%(s[0], str(int(s[1])+idxAdd), float(s[2]), float(s[3]))
        newline = formatted + "\n"
      elif "bond " in line:
        formatted = "%4s    %5s%5s   "%(s[0], str(int(s[1])+idxAdd), str(int(s[2])+idxAdd))
        newline = formatted + "   ".join(s[3:]) + "\n"
      elif "angle " in line:
        formatted = "%5s    %5s%5s%5s   "%(s[0], str(int(s[1])+idxAdd), str(int(s[2])+idxAdd), str(int(s[3])+idxAdd))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "anglep " in line:
        formatted = "%6s    %5s%5s%5s   "%(s[0], str(int(s[1])+idxAdd), str(int(s[2])+idxAdd), str(int(s[3])+idxAdd))
        newline = formatted + "   ".join(s[4:]) + "\n"
      elif "torsion " in line:
        formatted = "%7s    %5s%5s%5s%5s   "%(s[0], str(int(s[1])+idxAdd), str(int(s[2])+idxAdd), str(int(s[3])+idxAdd), str(int(s[4])+idxAdd))
        newline = formatted + "   ".join(s[5:]) + "\n"
      elif "multipole " in line:
        mpoleframes = s[1:-1]
        for i in range(len(mpoleframes)):
          if int(mpoleframes[i]) == 0:
            mpoleframes[i] = 0 
          elif int(mpoleframes[i]) < 0:
            mpoleframes[i] = "-"+str(abs(int(mpoleframes[i]))+idxAdd) 
          else:
            mpoleframes[i] = str(int(mpoleframes[i])+idxAdd) 
        if len(mpoleframes) == 4:# Z-Bisector or 3-Fold
          formatted = "%9s    %5s%5s%5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1], mpoleframes[2], mpoleframes[3])
        if len(mpoleframes) == 3:# Z-then-X or Bisector
          formatted = "%9s    %5s%5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1], mpoleframes[2])
        if len(mpoleframes) == 2:# Z-only
          formatted = "%9s    %5s%5s   "%(s[0], mpoleframes[0], mpoleframes[1])
        if len(mpoleframes) == 2:# Charge only 
          formatted = "%9s    %5s   "%(s[0], mpoleframes[0])
        newline = formatted + "         " + s[-1] + "\n"
      else:
        newline = line
      f.write(newline) 
  return

# Execute modKey to produce ttt.xyz and ttt.key
modKEY(key)
