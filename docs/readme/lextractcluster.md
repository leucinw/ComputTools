## lextractcluster.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lextractcluster.py)

* Usage:

	```shell
	python lextractcluster.py -i INPUT -r RADII [-n1 FIRST] [-n2 FINAL]
	```
* Intro:

	This program is used to extract molecular clusters from a big tinker xyz file, either a big cluster or a box. A `radii` is defined from the center of geometry of the solute molecule, which is defined by two atom numbers (`-n1` and `-n2`).
