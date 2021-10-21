# StarChart
**ManycoreOS Project** This is github repository for StarChart project from Seoul National University.

StarChart is a new graph-based similarity search algorithm for high-dimensional dataset. StarChart is an energy efficient algorithm that achieves lower query processing time and high throughput.

We will update the source code of StarChart when our paper is published.

StarChart has MIT license.

### Requirement List
1. Need dataset read wrapper for new dataset
2. Hardcoded \# of threads to utilize

### How to run
**Set up external library for optimized exection**

Clone external library from https://github.com/greg7mdp/parallel-hashmap which provides excellent hash map implementation. 
```bash
git clone https://github.com/greg7mdp/parallel-hashmap
```
**Set up StarChart**
```bash
git clone https://github.com/oslab-swrc/StarChart
cd StarChart
cp <path to parallel-hashmap>/parallel-hashmap ./src/lib -r
mkdir bin
make ssa
```

### Example output
Recall and total latency(ms) are reported as output

<img width="30%" src="example_output.png">
