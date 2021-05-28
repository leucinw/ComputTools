## ldimerinteract.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/ldimerinteract.py)

* Usage:

	```shell
	python ldimerinteract.py -m MOLECULES [MOLECULES ...] -p PROBES [PROBES ...]
	```
	
* Intro:

	This program is used to generate dimer structures for the given two monomers. In fact, it can support multiple monomers as input. For example, `python ldimerinteract.py -m Water.xyz Methanol.xyz -p Water.xyz Methanol.xyz`, will generate Water-Water, Water-Methanol, Methanol-Water (redundant) and Methanol-Methanol dimers. The atom types defined by MM2 are used to remove the degeneracy. This tend to give fewer structures. In the future, the molecular graph based symmetry will be applied.
