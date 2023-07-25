#!/bin/bash

# if [ $# -ne 1 ]; then
#     echo "Usage: $0 <mtx_file>"
#     exit 1
# fi
# mtx_file="$1"
# header=$(head -n 2 "$mtx_file")
# vertices=$2
# edges=$3



# source='1'
# insertedEdges='500'
# deletedEdges='0'
# inputGraph="Weighted$mtx_file"

# #convert unweighted mtx to random weighted mtx
# clang++ -std=c++17 weightedMtxConverter.cpp -o weightedMtxConverter && ./weightedMtxConverter -i $mtx_file -o $inputGraph -s 1 -e 10

# # Command 1 to generate random graph. v = # of vertices, e = # of edges, b = beginning of the range, t = terminating range, o = output file name
# clang++ -std=c++17 randomConnectedGraph.cpp -o randomConnectedGraph && ./randomConnectedGraph -v $vertices -e $edges -b 1 -t 10 -o $inputGraph

# # Command 2 to generate changedEdges. i for insertion count, d for deletion count, g for input graph name, o for output changed file.
# clang++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx

# # Command 3 to run sosp code. g for input graph, c for changed edges, s for source. 
# clang++ -std=c++17 SOSP.cpp -o program && ./program -g $inputGraph -c changeEdges.mtx -s $source

# # Command 4 to run parallelSOSP code, 
# #clang++ -std=c++11 -Xpreprocessor -fopenmp -lomp parallelSOSP.cpp -o parallelSOSP && ./parallelSOSP -g $inputGraph -c changeEdges.mtx -s $source







vertices=10000
edges=15000
source='1'
insertedEdges='15000'
deletedEdges='0'
inputGraph="weightedGraph.mtx"


# Command 1 to generate random graph. v = # of vertices, e = # of edges, b = beginning of the range, t = terminating range, o = output file name
#clang++ -std=c++17 randomConnectedGraph.cpp -o randomConnectedGraph && ./randomConnectedGraph -v $vertices -e $edges -b 1 -t 10 -o $inputGraph
cd /Users/smshovan/Documents/Codes/candy/PaRMAT-master/Release && make && ./PaRMAT -nVertices $vertices -nEdges $edges -sorted -noEdgeToSelf -noDuplicateEdges -threads 4 -output noWeightGraph.mtx && clang++ -std=c++17 noWeightedtoWeightedMtx.cpp  -o noWeightedtoWeightedMtx && ./noWeightedtoWeightedMtx -i noWeightGraph.mtx -o $inputGraph -s 1 -e 10 && cp $inputGraph /Users/smshovan/Documents/Codes/candy/tools/ && cd /Users/smshovan/Documents/Codes/candy/tools

# Command 2 to generate changedEdges. i for insertion count, d for deletion count, g for input graph name, o for output changed file.
clang++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx

# # Command 3 to run sosp code. g for input graph, c for changed edges, s for source. 
#clang++ -std=c++17 SOSP.cpp -o program && ./program -g $inputGraph -c changeEdges.mtx -s $source

# # Command 4 to run parallelSOSP code, 
# clang++ -std=c++11 -Xpreprocessor -fopenmp -lomp parallelSOSP.cpp -o parallelSOSP && ./parallelSOSP -g $inputGraph -c changeEdges.mtx -s $source

# # Command 5 to run parallelSOSPwithLoopSelection code parallelSOSPwithLoopSelection.cpp
# clang++ -std=c++11 -Xpreprocessor -fopenmp -lomp parallelSOSPwithLoopSelection.cpp -o parallelSOSPwithLoopSelection && ./parallelSOSPwithLoopSelection -g $inputGraph -c changeEdges.mtx -s $source

# # Command 6 to run parallelSOSPwithLoopSelectionWithoutCritical.cpp
 clang++ -std=c++11 -Xpreprocessor -fopenmp -lomp parallelSOSPwithLoopSelectionWithoutCritical.cpp -o parallelSOSPwithLoopSelectionWithoutCritical && ./parallelSOSPwithLoopSelectionWithoutCritical -g $inputGraph -c changeEdges.mtx -s $source
