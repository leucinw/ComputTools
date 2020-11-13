
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import os,sys,numpy,subprocess

usage = "python lGaussianJobRepair.py"

# color
RED = '\33[91m'
GREEN = '\33[92m'
YELLOW = '\33[93m'
ENDC = '\033[0m'

class Gaussian():
  def __init__(self, checkName):
    if ".log" in checkName:
      self.fname = checkName.split(".log")[0]
    elif ".com" in checkName:
      self.fname = checkName.split(".com")[0]
    else:
      sys.exit(RED + "provide a .com or .log file !" + ENDC)

    self.com = self.fname + ".com"
    self.log = self.fname + ".log"

  def getError(self):
    lines = open(self.log).readlines()
    err = None
    for line in lines:
      if "Error termination via" in line:
        err = line.split(".exe")[0].split("/")[-1]
        break
    return err
  
  def getCoords(self):
    babelcmd = "babel -ig09 %s -oxyz tmp.xyz && wait"%self.log
    subprocess.run(babelcmd, shell=True)
    xs,ys,zs = numpy.loadtxt("tmp.xyz", usecols=(1,2,3), skiprows=2, dtype="float", unpack=True)
    atoms = numpy.loadtxt("tmp.xyz", usecols=(0), skiprows=2, dtype="str", unpack=True)
    return atoms,xs,ys,zs

   
def main():
  files = os.listdir()
  for f in files:
    if f.endswith(".log"):
      obj = Gaussian(f)
      err = obj.getError()
      if err == "l9999": 
        [atoms, xs, ys, zs] = obj.getCoords()
        lines = open(obj.com).readlines()
        with open(obj.com + "_1", "w") as fcom:
          for i in range(0, len(lines)-len(xs)-1):
            if ("#" in lines[i]) and ("OPT" in lines[i].upper()):
              lines[i] = lines[i].upper().replace("OPT(", "OPT(CalcFC,").lower()
            fcom.write(lines[i])
          for i in range(0, len(xs)):
            fcom.write("%3s%12.6f%12.6f%12.6f\n"%(atoms[i], xs[i], ys[i],zs[i]))
          fcom.write("\n")
        os.rename(obj.com+"_1", obj.com)
        os.remove(obj.log)
        print(GREEN + "repaired %s with Error L9999 case"%f + ENDC)
      elif err == "l103": 
        [atoms, xs, ys, zs] = obj.getCoords()
        lines = open(obj.com).readlines()
        with open(obj.com + "_1", "w") as fcom:
          for i in range(0, len(lines)-len(xs)-1):
            if ("#" in lines[i]) and ("OPT" in lines[i].upper()):
              lines[i] = lines[i].upper().replace("OPT(", "OPT(Cartesian,").lower()
            fcom.write(lines[i])
          for i in range(0, len(xs)):
            fcom.write("%3s%12.6f%12.6f%12.6f\n"%(atoms[i], xs[i], ys[i],zs[i]))
          fcom.write("\n")
        os.rename(obj.com+"_1", obj.com)
        os.remove(obj.log)
        print(GREEN + "repaired %s with Error L103 case"%f + ENDC)
      elif err == "l202": 
        [atoms, xs, ys, zs] = obj.getCoords()
        lines = open(obj.com).readlines()
        with open(obj.com + "_1", "w") as fcom:
          for i in range(0, len(lines)-len(xs)-1):
            if ("#" in lines[i]) and ("Symmetry=None" not in lines[i]): 
              lines[i] = lines[i][:-1] + " Symmetry=None \n"
            fcom.write(lines[i])
          for i in range(0, len(xs)):
            fcom.write("%3s%12.6f%12.6f%12.6f\n"%(atoms[i], xs[i], ys[i],zs[i]))
          fcom.write("\n")
        os.rename(obj.com+"_1", obj.com)
        os.remove(obj.log)
        print(GREEN + "repaired %s with Error L202 case"%f + ENDC)
      elif err == None:
        pass
      else:
        print(RED + "Can not repair %s currently"%err + ENDC)

if __name__ == "__main__":
  main() 
    
