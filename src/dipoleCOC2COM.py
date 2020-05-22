
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import numpy as np
import sys

gaufile = sys.argv[1]

lines = open(gaufile).readlines()

for line in lines:
  if "Standard orientation" in line:
    idx = lines.index(line) + 5
  if "Sum of Mulliken charges" in line:
    charge = float(line.split()[-1])
  if " X= " in line:
    data = line.split()
    dpx = float(data[1]) 
    dpy = float(data[3])
    dpz = float(data[5])
    tot = float(data[7])

x=[];y=[]; z=[]; element=[]
while True:
  if "----" in lines[idx]:break
  data = lines[idx].split()
  x.append(float(data[-3]))
  y.append(float(data[-2]))
  z.append(float(data[-1]))
  element.append(data[1])
  idx += 1

x = np.array(x)
y = np.array(y)
z = np.array(z)

nucDict = {"6":6.0 , "1":1.0, "8":8.0, "7":7.0, "15":15.0, "16":16.0}
massDict = {"6":12.011 , "1":1.00797, "8":15.9994, "7":14.0067, "15":30.973762, "16":32.0}

#COM
com_x = 0.0
com_y = 0.0
com_z = 0.0
masssum = 0.0
nucsum = 0.0
for i in range(len(element)):
  masssum += massDict[element[i]] 
  nucsum += nucDict[element[i]]

for i in range(len(x)):
  com_x = com_x + x[i]*massDict[element[i]]/masssum 
  com_y = com_y + y[i]*massDict[element[i]]/masssum
  com_z = com_z + z[i]*massDict[element[i]]/masssum

COM = np.array([com_x,com_y,com_z] )

#COC
coc_x = 0.0
coc_y = 0.0
coc_z = 0.0
for i in range(len(x)):
  coc_x = coc_x + x[i]*nucDict[element[i]]/nucsum 
  coc_y = coc_y + y[i]*nucDict[element[i]]/nucsum
  coc_z = coc_z + z[i]*nucDict[element[i]]/nucsum
COC = np.array([coc_x,coc_y,coc_z])

bohr_e2debye = 0.393456 
bohr2A = 0.529177210 
vec_gaus = charge*(COC-COM)/bohr2A/bohr_e2debye
dipole_COC = np.array([dpx, dpy, dpz])

dipole_COM = dipole_COC + vec_gaus
dipole_COM_t = ((dipole_COM)**2).sum()**0.5
print("%10s%10s%10s%10s%10s"%(" ", "X", "Y", "Z", "Tot"))
print("%10s%10.4f%10.4f%10.4f%10.4f"%("COC",dipole_COC[0], dipole_COC[1], dipole_COC[2], tot))
print("%10s%10.4f%10.4f%10.4f%10.4f"%("COM",dipole_COM[0], dipole_COM[1], dipole_COM[2], dipole_COM_t))
