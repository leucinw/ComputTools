
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import numpy as np
from scipy.optimize import least_squares 

def f_diff(params):
  QM, weight = np.loadtxt('QM-energy.dat',usecols=(1,2))
  MM = np.array(getMM(params))
  indicator = np.square(QM-MM).mean()
  print(" Current MSE is %15.8f !"%indicator, end="\r")
  return (QM-MM)*weight

def main():
  x0 = np.loadtxt("p0.txt",usecols=(-1,), dtype='float')
  myBounds = ([0.3]*len(x0), [2.5]*len(x0)) # use bounds
  ret = least_squares(f_diff, x0, bounds=myBounds)
  np.savetxt("p1.txt", ret.x,fmt='%15.10f')

if __name__ == "__main()__":
  main()
