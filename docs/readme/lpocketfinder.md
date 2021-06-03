## lpocketfinder.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lpocketfinder.py)

* Usage:

	```shell
	python lpocketfinder.py protein.pdb ligand.mol2
	```

* Intro:

	This program is used to find the possible binding pocket of a protein using a ligand as the probe. This is a brute-force scan of the 3D grid of a protein structure using GOLD docking software. Besides the protein.pdb and ligand.mols files, two additional files for running GOLD are required:

	* gold.conf: regular gold.conf file 
	* sub.sh: should be $GOLD_PATH/gold_auto gold.conf
