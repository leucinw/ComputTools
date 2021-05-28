## lsub.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/lsub.py)

* Usage:

	```shell
	python lsub.py -i INPUT [INPUT ...] [-n NODES [NODES ...]] [-x XNODES [XNODES ...]] [-d DISK] [-m MEMORY] [-M MAXMEM] [-t TCHECK]
	```

* Intro:

	This program is used to submit QM jobs onto renlab clusters. Currently it supports gaussian (.com), psi4 (.psi4) and qchem (.qchem) formats. It is very flexible to specify/exclude the required computer nodes. `ssh` is used to spread the jobs to CPU nodes. The program will not terminate until all specified jobs have been submitted.

* Examples:

	It allows very flexible node selection. Here are several examples:

	1. Submit all .com files; check the node availability every 30 second; on >50 GB memory nodes

	```shell
	python lsub.py -i *.com -m 50 -t 30
	```

	1. Submit all .com files; check the node availability every 30 second; on node144, node142 only 

	```shell
	python lsub.py -i *.com -n node142 node144 -t 30
	```
	1. Submit all .com files; check the node availability every 30 second; on >50 GB memory nodes but exclude node142 and node144

	```shell
	python lsub.py -i *.com -x node142 node144 -t 30 -m 50
	```
