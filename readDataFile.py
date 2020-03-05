
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


def readWholeFile(filename):
    data = []
    for line in open(filename).readlines():
      if len(line) != 0:
        data.append(line.split("\n")[0])
    return data

def readOneColumn(filename,ncolum):
    data=[]
    for line in open(filename).readlines():
      strings=line.split()[ncolum]
      data.append(strings)
    return data

def readTwoColumns(filename):
    xdata=[]
    ydata=[]
    for line in open(filename).readlines():
      xdata.append(float(line.split()[0]))
      ydata.append(float(line.split()[1]))
    return xdata,ydata

def readThreeColumns(filename):
    xdata=[]
    ydata=[]
    zdata=[]
    for line in open(filename).readlines():
      xdata.append(float(line.split()[0]))
      ydata.append(float(line.split()[1]))
      zdata.append(float(line.split()[2]))
    return xdata,ydata,zdata

def readFourColumns(filename):
    xdata=[]
    ydata=[]
    zdata=[]
    wdata=[]
    for line in open(filename).readlines():
      xdata.append(float(line.split()[0]))
      ydata.append(float(line.split()[1]))
      zdata.append(float(line.split()[2]))
      wdata.append(float(line.split()[3]))
    return xdata,ydata,zdata,wdata

def readKeywordLines(filename,keyword):
    xdata=[]
    for line in open(filename).readlines():
      if keyword in line:
        xdata.append(line)
    return xdata
