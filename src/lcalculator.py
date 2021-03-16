
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import argparse
import numpy as np

# color
RED = '\033[91m'
GREEN = '\033[92m'
ENDC = '\033[0m'

def main():
  #===>>>
  parser = argparse.ArgumentParser()
  parser.add_argument('-f1', dest = 'datafile1', required=True)  
  parser.add_argument('-c1', dest = 'ncolumn1', required=True, type=int)  
  parser.add_argument('-f2', dest = 'datafile2', default=None)  
  parser.add_argument('-c2', dest = 'ncolumn2',required=False, type=int)  
  parser.add_argument('-s1', dest = 'skiprow1', default=0, type=int)  
  parser.add_argument('-s2', dest = 'skiprow2', default=0, type=int)  
  parser.add_argument('-o',  dest = 'operator', required=True, type=str.upper, choices=['MIN', 'MAX', 'MEDIAN', 'MEAN', 'STD', 'RMSE', 'MUE', 'MSE'])
  args = vars(parser.parse_args())
  
  #===>>> 
  f1 = args["datafile1"]
  f2 = args["datafile2"]
  c1 = args["ncolumn1"]
  if args["ncolumn2"]:
    c2 = args["ncolumn2"]
  op = args["operator"]
  s1 = args["skiprow1"]
  s2 = args["skiprow1"]
  
  # minimal
  def fmin(arr):
    print(GREEN + "MIN value: %12.6f"%(np.min(arr)) + ENDC)
    return 
  # maximal
  def fmax(arr):
    print(GREEN + "MAX value: %12.6f"%(np.max(arr)) + ENDC)
    return 
  # median
  def fmedian(arr):
    print(GREEN + "MEDIAN value: %12.6f"%(np.median(arr)) + ENDC)
    return 
  # mean
  def fmean(arr):
    print(GREEN + "MEAN value: %12.6f"%(np.mean(arr)) + ENDC)
    return 
  # standard deviation
  def fstd(arr):
    print(GREEN + "Standard Deviation: %12.6f"%(np.std(arr)) + ENDC)
    return 
  # root mean square error/deviation
  def frmse(arr1,arr2):
    rmse = np.sqrt(np.square(arr1-arr2).mean())
    print(GREEN + "Root Mean Square Error: %12.6f"%(rmse) + ENDC)
    return 
  # mean signed error/deviation
  def fmse(arr1,arr2):
    mse = (arr1-arr2).mean()
    print(GREEN + "Mean Signed Error: %12.6f"%(mse) + ENDC)
    return 
  # mean unsigned error/deviation
  def fmue(arr1,arr2):
    mue = (np.abs(arr1-arr2)).mean()
    print(GREEN + "Mean Unsigned Error: %12.6f"%(mue) + ENDC)
    return 
  
  arr1 = np.loadtxt(f1, usecols=(c1,), dtype="float", skiprows=s1)
  if f2:
    arr2 = np.loadtxt(f2, usecols=(c2,), dtype="float", skiprows=s2)

  if (op == "MIN"):
    fmin(arr1)
  elif (op == "MAX"):
    fmax(arr1)
  elif (op == "MEAN"):
    fmean(arr1)
  elif (op == "MEDIAN"):
    fmedian(arr1)
  elif (op == "STD"):
    fstd(arr1)
  elif (op == "RMSE"):
    frmse(arr1, arr2)
  elif (op == "MSE"):
    fmse(arr1, arr2)
  elif (op == "MUE"):
    fmue(arr1, arr2)
  else:
    print(RED + "Operator '%s' not supported!"%op + ENDC)
  return

if __name__ == "__main__":
  main()
