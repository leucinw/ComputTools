## ltorsion.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/ltorsion.py)

* Usage:

  ```shell
  python ltorsion.py -xyz INPUTXYZ -key INPUTKEY [-optm OPTMETHOD]
	                   [-optb OPTBASIS] [-spm SPMETHOD] [-spb SPBASIS]
										 [-chrg CHARGE] [-spin SPIN] -mode MODE [-npt NPOINT]
										 [-intv INTERVAL] [-rstr RESTRAIN] [-rig RIGIDFIT]
  ```

* Intro:

  This program is used to (1) generate QM and MM files for torsion scan; (2) fitting torsional parameters for MM model. This is part of `prmGen` program but can be used alone.