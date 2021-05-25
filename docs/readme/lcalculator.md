## lcalculator.py

* Source: 
	
	[https://github.com/leucinw/ComputTools/tree/master/src/lcalculator.py](https://github.com/leucinw/ComputTools/tree/master/src/lcalculator.py)

* Usage:

	```shell
	python lcalculator.py -f1 DATAFILE1   -c1 NCOLUMN1 
											 [-f2 DATAFILE2] [-c2 NCOLUMN2] 
											 [-s1 SKIPROW1]  [-s2 SKIPROW2] 
											 -o {MIN,MAX,MEDIAN,MEAN,STD,SUM,RMSE,MUE,MSE}
	```

* Intro:
	
	This is a calculator that can deals with statistics such as min, max, median, RMSE etc. Options are explained below:
	
	* -f1: first data file [required]
	* -f2: second data file [optional]
	* -c1: column number of the first file
	* -c2: colunm number of the second file
	* -s1: skip the first `s1` rows 
	* -s2: skip the first `s2` rows
	* -o: operation to take
