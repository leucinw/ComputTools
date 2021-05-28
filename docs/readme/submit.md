## submit.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/submit.py)

* Usage:

	```shell
	python submit.py [-p PATH] -x SCRIPTS [SCRIPTS ...]
	```
* Intro:

	This program is used to submit one or multiple script(s) to one CPU node. One must login to a computing node to execute this program. This is very useful for submitting light jobs such as `analyze`, `minimize`, `polarize` etc of Tinker software. `concurrent` module of Python is used to parallellize the jobs.

* Example:
	
	* submit files in the current folder
	```shell
	python submit.py -x ana1.sh ana2.sh ana3.sh ... ana1000.sh 
	```

	* submit files in another folder
	```shell
	python submit.py -p /home/liuchw/ -x ana1.sh ana2.sh ana3.sh ... ana1000.sh 
	```

