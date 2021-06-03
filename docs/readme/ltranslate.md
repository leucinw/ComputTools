## ltranslate.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/ltranslate.py)

* Usage:

  ```shell
  python ltranslate.py -xyz XYZFILE -p1 POINTS1 [POINTS1 ...] -p2 POINTS2 [POINTS2 ...] -n1 NFIRST [-d DISTANCE]
  ```

* Intro:

  This program is used to generate a series of dimers with different intermolecular distances based on the given geometry. The translation vector can be defined by (1) two points (one on each monomer), (2) geometry centers of the atoms defined in each monomer. The second translation vector definition enables us to translate special configurations, such as pi-pi stacking. 
