#!/usr/bin/env python

import string

def readWholeFile(filename):
    data = []
    lines = open(filename).readlines()
    for line in lines:
      if len(line) != 0:
        data.append(line.split("\n")[0])
    return data

def readOneColumn(filename,ncolum):
    data=[]
    lines=open(filename).readlines()
    for n in range(len(lines)):
      strings=lines[n].split()[ncolum]
      if "D" in strings:
        strings=strings.replace("D", "e")
      data.append(float(strings))
    return data

def readTwoColumns(filename):
    xdata=[]
    ydata=[]
    lines=open(filename).readlines()
    for n in range(len(lines)):
      xdata.append(float(lines[n].split()[0]))
      ydata.append(float(lines[n].split()[1]))
    return xdata,ydata

def readThreeColumns(filename):
    xdata=[]
    ydata=[]
    zdata=[]
    lines=open(filename).readlines()
    for n in range(len(lines)):
      xdata.append(float(lines[n].split()[0]))
      ydata.append(float(lines[n].split()[1]))
      zdata.append(float(lines[n].split()[2]))
    return xdata,ydata,zdata

def readFourColumns(filename):
    xdata=[]
    ydata=[]
    zdata=[]
    wdata=[]
    lines=open(filename).readlines()
    for n in range(len(lines)):
      xdata.append(float(lines[n].split()[0]))
      ydata.append(float(lines[n].split()[1]))
      zdata.append(float(lines[n].split()[2]))
      wdata.append(float(lines[n].split()[3]))
    return xdata,ydata,zdata,wdata

def readKeywordLines(filename,keyword):
    xdata=[]
    lines=open(filename).readlines()
    for line in lines:
      if keyword in line:
        xdata.append(line)
    return xdata

