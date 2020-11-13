

# -E(l) + E(g) + RT
import numpy as np
E_liquid = []
E_gas = []

for line in open("liquid.log").readlines():
  if "Total Potential Energy" in line:
    E_liquid.append(float(line.split()[-2]))

for line in open("gas.log").readlines():
  if "Total Potential Energy" in line:
    E_gas.append(float(line.split()[-2]))

natom_liquid = int(open("liquid.xyz").readlines()[0].split()[0])
natom_gas = int(open("gas.xyz").readlines()[0].split()[0])

nMol = natom_liquid/natom_gas
R = 1.98720425864083 * 0.001 
T = 174.15
kcal2kjol = 4.184
Hvap = (- np.array(E_liquid).mean()/nMol + np.array(E_gas).mean() + R*T)*kcal2kjol

print(Hvap)
