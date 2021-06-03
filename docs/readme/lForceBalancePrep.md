## lForceBalancePrep.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lForceBalancePrep.py)

* Usage:

	```shell
	python lForceBalancePrep.py input.dat
	```

* Intro:

	This program is used to prepare an `interaction.txt` that is used as a ForceBalance interaction energy target. A typical `input.dat` file contains three columns
	
	```
	Filename       MP2/CBS       Weight
	xxx_0.70.xyz   -5.00         0.5
	xxx_0.80.xyz   -4.00         0.6
	xxx_0.90.xyz   -3.00         0.7
	xxx_0.95.xyz   -2.00         0.8
	xxx_1.00.xyz   -1.00         1.0
	```
	
	The calculation of dimer interaction energy in ForceBalance requires dimer and both monomers. A program called [lSplitCluster.py](./lSplitCluster.md) is used to split a dimer into two monomers.
