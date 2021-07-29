## lChargeIntegerizer.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/lChargeIntegerizer.py)

* Usage:

  ```shell
  python lChargeIntegerizer.py -xyz XYZ -prm PRM -chg CHG
  ```

* Intro:

  This program is used to make the net charge of the molecule an integer. Sometimes a molecular fragment may have very small amount of non-integer charge, in which case we need to make it an integer in order to do molecular dynamics simulations. The approach used here is to distribute the small amount of residual charge onto the heavy atoms of the molecule, with the intention of avoiding overcharging of some atoms.