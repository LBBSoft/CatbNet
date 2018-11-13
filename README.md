# CatbNet

Tools for comparing and analyzing the topology of biological networks

## Getting Started
CatbNet is able to compute most important network topological features, from node-based measures such as: different centralities to whole network-based measures such as: network density, diameter, for analytical and comparative purposes. 
You can use this tool by downloading its python source code.
This tool is created using python 2.7

### Prerequisites

If you are interested in using python source code, consider installing packages:
```
Networkx  version 1.9.1
Numpy
Matplotlib
Pandas
Scipy
PyQt4
```
## Running and Loading Data
To run CatbNet from python source run command:

```
python main.py
```

In CatbNet there is the opportunity to laod data from .gml, .graphml, .net and .el network data formats.
In the case of loading .el (edge list) networks you should consider that:
* Each line must demonstrate a node or an edge
* Data must be tab seperated
* In weighted networks, each edge line must have three values in a line
* Nodes and Edges are seprated by #nodes and #edges keywords

```
example:

#nodes
a
b
c
#edges
a	b	1.45
a	c	2.6
```

### Considerations

In the case you want to compare network groups, you should check option 'network files are classified in groups'.


### Test Application

As a Test Data, you can test network files of grouped networks in directory: test-data


## License

This project is licensed under the LBB License - It is free to use for only researching goals.
For more details contact:

```
Ehsan Pournoor

Email: e.pournoor@ut.ac.ir
Laboratory of Bioinformatics and Systems Biology (LBB)
```

