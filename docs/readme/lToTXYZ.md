## lToTXYZ.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/lToTXYZ.py)

* Usage:

  ```shell
  python lToTXYZ.py -i INP -t {xyz,g09,qcout,mol,mol2,psi4,sdf,pdb,psi4out}
  ```

* Intro:

  This program is used to convert other files to tinker xyz format. Note that it is different from `babel` provided functionality, which uses MM2 style atom type. In this program, I use molecular graph based symmetry to define the atom types. 