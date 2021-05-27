## simBARsim.py

* Source: 

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/simBARsim.py)

* Usage:
	
	```shell
	python simBARsim.py option(s)
	```

* Intro:

	The free energy simulations using BAR program involve multiple windows and steps. Here I provide this script to simplify this process. This program runs with 4 options, either individually or together: 

	* setup: generate windows and necessary files
	* dynamic: parallelly run molecular dynamics simulation with either CPU or GPU (set in dynamic.sh)
	* bar: parallelly run bar simulation to get .bar and free energy (set in bar.sh)
	* result: get the final result and write in `result.txt` file

	In practical simulations, the following command is allowed:

	```shell
	python simBARsim.py setup dynamic bar result
	```
	
* Example:

	An example folder to run K+ solvation in water is provided [here](https://github.com/leucinw/ComputTools/tree/master/docs/data/bardemo)
