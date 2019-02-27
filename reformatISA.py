#!/usr/bin/env python

# Usage: reformatISA.py isa_multipoles.txt

# Convert the multipole from Iterative Stockholder Atoms (ISA)
# to the GDMA-like format in order to be read by TINKER poledit

# Author: Chengwen Liu
# University of Texas at Austin
# June 9th, 2017

import os,sys

#Standard ISA output file from HiPart
isaOUT = sys.argv[1]

#Gaussian output file for Cartesian
cartesian = sys.argv[2]

#Formatted ISA multipole file for use with TINKER
formattedISA = sys.argv[3]

oFile = open(formattedISA,'w')
strings = 20*" "
oFile.write(strings+"Iterative Stockholder Atoms\n\n")
oFile.write("ISA from HiPart\n\n")
oFile.write("Positions and radii in angstrom\n")
oFile.write("Multipole moments in atomic units, ea_0^k for rank k\n\n")

#Read ATOM and COORDINATES from Gaussian 
atomList = []
coordList = []
lines = file(cartesian).readlines()
for line in lines:
  if "Standard orientation" in line:
    nIndex = lines.index(line)
i = 5
while True:
  data = lines[nIndex+i].split()
  if len(data)==6: 
    coordList.append([float(data[3]), float(data[4]), float(data[5])])
    i+=1
  else:
    break

#Read Multipoles from ISA output
multipole = []
lines = file(isaOUT).readlines()
for line in lines:
  if ( "|" in line) and ("Multipoles" not in line):
    rawData = line.split()
    atomList.append(rawData[1])
    Q00 = float(rawData[4])
    Q10 = float(rawData[5]); Q11c=float(rawData[6]); Q11s=float(rawData[7])
    Q20 = float(rawData[8]); Q21c=float(rawData[9]); Q21s=float(rawData[10]); Q22c = float(rawData[11]); Q22s=float(rawData[12])
    multipole.append([Q00, Q10, Q11c, Q11s, Q20, Q21c, Q21s, Q22c, Q22s])

#Write them to the formatted multipole file

for n in xrange(len(atomList)):
  oFile.write("%-11sx =%10.6f  y =%10.6f  z =%10.6f\n\n"%(atomList[n], coordList[n][0], coordList[n][1], coordList[n][2]))
  oFile.write("%22s  = %10.6f\n"%("Q00", multipole[n][0]))
  oFile.write("%22s  = %10.6f%6s =%10.6f%7s =%10.6f\n"%("Q10", multipole[n][1], "Q11c", multipole[n][2], "Q11s", multipole[n][3]))
  oFile.write("%22s  = %10.6f%6s =%10.6f%7s =%10.6f\n"%("Q20", multipole[n][4], "Q21c", multipole[n][5], "Q21s", multipole[n][6]))
  oFile.write("%23s = %10.6f%6s =%10.6f\n"%("Q22c", multipole[n][7], "Q22s", multipole[n][8],))
oFile.close()

