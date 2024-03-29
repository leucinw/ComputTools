## subsubmit.py

* Source:

	[click here](https://github.com/leucinw/ComputTools/tree/master/src/subsubmit.py)

* Usage:

	```shell
	python subsubmit.py -x SCRIPTS [SCRIPTS ...] [-n NODES [NODES ...]]
	```
* Intro:

	This program is used to submit one or multiple script(s) to one or multiple CPU node. One does not have to login to a computing node to execute this program. This is very useful for submitting light jobs such as `analyze`, `minimize`, `polarize` etc of Tinker software. `parallel` tool of GNU is used to submit the high level script files. On a certain node, `concurrent.futures` model of Python is used to submit the low-level scripts files, as describe in [submit](./submit.md).

* Example:
	
	```shell
	python subsubmit.py -x ana1.sh ana2.sh ... ana1000.sh -n node90 node91 ... node100
	```

* Efficiency:
	
	Here is the efficiency validations on 2163 `polarize` jobs.

	| Command    | Description | Time(s)|
	| ----------- | ----------- |--------|
	| `sh polar*.sh`      | sequencial       |   396     |
	| `subsubmit.py -x polar*.sh -n node91`  | parallel; 1 node        | 36|
  | `subsubmit.py -x polar*.sh -n node91 node...`| parallel; 2 nodes | 20|
  | `subsubmit.py -x polar*.sh -n node91 node...`| parallel; 4 nodes | 11|
  | `subsubmit.py -x polar*.sh -n node91 node...`| parallel; 6 nodes | 8 |
  | `subsubmit.py -x polar*.sh -n node91 node...`| parallel; 8 nodes | 7 |

	As seen above, more than 50x accelaration is obtained by using both `GNU parallel` and `concurrent`. Depending on how many CPU cores your `submitting` node has, the efficiency may vary. Consulting [GNU parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) for more information.
	

