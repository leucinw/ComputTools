
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

def readPsi4Out(OUT):
  for line in open(OUT).readlines():
  ele  = 0.0
  exch = 0.0
  ind  = 0.0
  disp = 0.0
  tot  = 0.0
  for line in lines:
    if "Electrostatics     " in line:
      ele = float(line.split()[3])
    elif "Exchange     " in line:
      exch = float(line.split()[3])
    elif "Induction     " in line:
      ind = float(line.split()[3])
    elif "Dispersion     " in line:
      disp = float(line.split()[3])
  vdw = exch + disp
  tot = vdw + ind + ele
  return ele, ind, exch, disp, vdw, tot
  
def readTinkerOut(OUT):
  ele  = 0.0
  vdw = 0.0
  pol = 0.0
  tot  = 0.0
  for line in open(OUT).readlines():
    if "Electrostatics     " in line:
      ele = float(line.split()[3])
    elif "Exchange     " in line:
      exch = float(line.split()[3])
    elif "Induction     " in line:
      ind = float(line.split()[3])
    elif "Dispersion     " in line:
      disp = float(line.split()[3])
  vdw = exch + disp
  tot = vdw + ind + ele
  return ele, ind, exch, disp, vdw, tot

