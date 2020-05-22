
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================


import argparse
import numpy as np

def main():
  #===>>>
  parser = argparse.ArgumentParser()
  parser.add_argument('-f1', dest = 'datafile1', required=True)  
  parser.add_argument('-f2', dest = 'datafile2', default=None)  
  parser.add_argument('-c1', dest = 'ncolumn1', required=True)  
  parser.add_argument('-c2', dest = 'ncolumn2')  
  parser.add_argument('-o',  dest = 'operator', required=True, choices=['MIN', 'MAX', 'MEDIAN', 'MEAN', 'STD', 'RMSE', 'MUE', 'MSE'])
  args = vars(parser.parse_args())
  
  #===>>> 
  f1 = args["datafile1"]
  f2 = args["datafile2"]
  c1 = int(args["ncolumn1"])
  c2 = int(args["ncolumn2"])
  op = args["operator"]
  
  # minimal
  def fmin(arr):
    print("MIN value: %12.6f"%(np.min(arr)))
    return 
  # maximal
  def fmax(arr):
    print("MAX value: %12.6f"%(np.max(arr)))
    return 
  # median
  def fmedian(arr):
    print("MEDIAN value: %12.6f"%(np.median(arr)))
    return 
  # mean
  def fmean(arr):
    print("MEAN value: %12.6f"%(np.mean(arr)))
    return 
  # standard deviation
  def fstd(arr):
    print("Standard Deviation: %12.6f"%(np.std(arr)))
    return 
  # root mean square error/deviation
  def frmse(arr1,arr2):
    rmse = np.sqrt(np.square(arr1-arr2).mean())
    print("Root Mean Square Error: %12.6f"%(rmse))
    return 
  # mean signed error/deviation
  def fmse(arr1,arr2):
    mse = (arr1-arr2).mean()
    print("Mean Signed Error: %12.6f"%(mse))
    return 
  # mean unsigned error/deviation
  def fmue(arr1,arr2):
    mue = (np.abs(arr1-arr2)).mean()
    print("Mean Unsigned Error: %12.6f"%(mue))
    return 
  
  arr1 = np.loadtxt(f1, usecols=(c1,), dtype="float", skiprows=0)
  if f2:
    arr2 = np.loadtxt(f2, usecols=(c2,), dtype="float", skiprows=0)

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
    print("Operator '%s' not supported!"%op)
  return

if __name__ == "__main__":
  main()
