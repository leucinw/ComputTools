## lAtomTypeGenerator.py

* Source:

  [click here](https://github.com/leucinw/ComputTools/tree/master/src/lAtomTypeGenerator.py)

* Usage:

  ```shell
  python lAtomTypeGenerator.py -i txyz -m mode
  ```

* Intro:

  This program is used to generate atom string types used in parameter assignment. These string types are `graph-based` and should be more robust comparing to the currently used SMARTS pattern match. Similar idea has been implemented in the paper [by Savoie and co-workers](https://chemrxiv.org/articles/preprint/Topology_Automated_Force-Field_Interactions_TAFFI_A_Framework_for_Developing_Transferable_Force-Fields/14527299). It support (1) mode 0, with 0 level of graph, which equals to element based type and (2) mode 1, with 1 level of neighboring atoms. Hydrogen is further split with addition of the number of connections of its heavy atom. For example, the atom types generated for Ethanol molecules will be:

	| atom No. | atom | atom type | string type |
	| ------|---  | -----     | -----       |
	| 1 | O | 12  |  OCH      | 
	| 2 | C | 10  |  CCHHO    | 
	| 3 | C | 11  |  CCHHH    | 
	| 4 | H | 7   |  HC4      | 
	| 5 | H | 7   |  HC4      | 
	| 6 | H | 8   |  HC4      | 
	| 7 | H | 8   |  HC4      | 
	| 8 | H | 8   |  HC4      | 
	| 9 | H | 9   |  HO2      | 
