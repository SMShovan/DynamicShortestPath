#!/bin/bash

source='1'
vertices='50000'
edges='700000'
insertedEdges='50000'
deletedEdges='0'
inputGraph='randomgraph.mtx'

# Command 1 to generate random graph. v = # of vertices, e = # of edges, b = beginning of the range, t = terminating range, o = output file name
clang++ -std=c++11 randomConnectedGraph.cpp -o randomConnectedGraph && ./randomConnectedGraph -v $vertices -e $edges -b 1 -t 10 -o $inputGraph

# Command 2 to generate changedEdges. i for insertion count, d for deletion count, g for input graph name, o for output changed file.
clang++ -std=c++11 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx

# Command 3 to run sosp code. g for input graph, c for changed edges, s for source. 
clang++ -std=c++11 SOSP.cpp -o program && ./program -g $inputGraph -c changeEdges.mtx -s $source
