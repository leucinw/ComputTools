## lChemFileEditor.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/lChemFileEditor.py)

* Usage:

  ```shell
  python lChemFileEditor.py -i INPUT -m MODE 
  ```

* Intro:

  This program is used to edit important computational chemistry files. Currently it supports tinker xyz and PDB file. 

	* For a tinker xyz file, it moves the Hydrogens to the attached heavy atoms (nicer mode) or splits a cluster file into monomers (split mode). 
	* For a PDB file, it splits the PDB file into its residues (split mode) or re-write the residues into a nicer PDB file (nicer mode), which will have atom number and residue number in good order.