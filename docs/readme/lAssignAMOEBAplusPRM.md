## lAssignAMOEBAplusPRM.py

* Source:

 	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lAssignAMOEBAplusPRM.py)

* Usage:

	```shell
	python lAssignAMOEBAplusPRM.py -i XYZ -k KEY -p {POLAR,CF,BONDED,NONBONDED}
	```

* Intro:
	
	This program is used to assign AMOEBA+ parameters for a new molecule. These parameters including `atomic polarizability`, `charge flux`, `bonded parameters` and `nonbonded parameters`. Currently, it will only assign the parameters if the chemistry is found in our existing database.
