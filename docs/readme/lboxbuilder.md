## lboxbuilder.py

* Source: 

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lboxbuilder.py)

* Usage:

	```shell
	python lboxbuilder.py -mode MODE -solute SOLUTE [-solvent SOLVENT] -prm PRM -size SIZE -density DENSITY -tinker PATH [-molar1 MOLAR1]
	```

* Intro:

	This program is used to generate simulation box for Tinker/Tinker-OpenMM/Tinker9. Currently the following functions are provided: 
	
	* Pure liquid box with specified density/box size
	* One molecule in water (for hydration) or other solvents (solvation) with specified density/box size
	* Binary mixture (need to provide the molar fraction of the first molecule) with specified density/box size

* Examples:
	* Generate a pure liquid box
	```shell
	python lboxbuilder.py -mode 1 -solute water.xyz -prm water03.prm -size 30 -density 0.997 -tinker $TINKER8
	```

	* Generate a methanol-in-water box
	```shell
	python lboxbuilder.py -mode 2 -solute methanol.xyz -solvent water.xyz -prm amoeba09.prm -size 30 -density 0.998 -tinker $TINKER8
	```

	* Generate a binary mixture of 50% methanol and 50% water
	```shell
	python lboxbuilder.py -mode 3 -solute methanol.xyz -solvent water.xyz -prm amoeba09.prm -size 30 -density 0.980 -tinker $TINKER8 -molar1 0.5
	```

