## matchPRM.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/matchPRM.py)

* Usage:

  ```shell
	python matchPRM.py template.prm template.txyz dealwith.txyz energyterm
  ```

* Intro:

  This program is used to assign the force field (AMOEBA/AMOEBA+) parameters for the `dealwith.txyz` according to the `template.txyz` and `template.prm`. The template and dealwith files may have different atom types or even different atom orders.