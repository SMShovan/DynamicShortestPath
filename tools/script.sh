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






#For random graph
vertices=10
edges=15


source='1'
insertedEdges='5'
deletedEdges='0'
inputGraph="weightedGraph.mtx"
dataset="inf-roadNet-CA.mtx"
threads="4"

## Real graph

#clang++ -std=c++17  mtxDatasetToRandomWeight.cpp -o mtxDatasetToRandomWeight && ./mtxDatasetToRandomWeight -i rgg_n_2_20_s0.mtx -o output_file_name.mtx -s 1 -e 10

# # Command 1 to generate random graph. v = # of vertices, e = # of edges, b = beginning of the range, t = terminating range, o = output file name
# #clang++ -std=c++17 randomConnectedGraph.cpp -o randomConnectedGraph && ./randomConnectedGraph -v $vertices -e $edges -b 1 -t 10 -o $inputGraph
# cd /Users/smshovan/Documents/Codes/candy/PaRMAT-master/Release && make && ./PaRMAT -nVertices $vertices -nEdges $edges -sorted -noEdgeToSelf -noDuplicateEdges -threads 4 -output noWeightGraph.mtx && clang++ -std=c++17 noWeightedtoWeightedMtx.cpp  -o noWeightedtoWeightedMtx && ./noWeightedtoWeightedMtx -i noWeightGraph.mtx -o $inputGraph -s 1 -e 10 && cp $inputGraph /Users/smshovan/Documents/Codes/candy/tools/ && cd /Users/smshovan/Documents/Codes/candy/tools

# # Command 2 to generate changedEdges. i for insertion count, d for deletion count, g for input graph name, o for output changed file.
# clang++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx

# # # Command 3 to run sosp code. g for input graph, c for changed edges, s for source. 
# clang++ -std=c++17 SOSP.cpp -o program && ./program -g $inputGraph -c changeEdges.mtx -s $source

# # # Command 6 to run parallelSOSPwithLoopSelectionWithoutCritical.cpp
# clang++ -std=c++11 -Xpreprocessor -fopenmp -lomp parallelSOSPwithLoopSelectionWithoutCritical.cpp -o parallelSOSPwithLoopSelectionWithoutCritical && ./parallelSOSPwithLoopSelectionWithoutCritical -g $inputGraph -c changeEdges.mtx -s $source


#Froundy

# #Real graph
#g++ -std=c++17  mtxDatasetToRandomWeight.cpp -o mtxDatasetToRandomWeight && ./mtxDatasetToRandomWeight -i $dataset -o $inputGraph -s 1 -e 10

#generate random graph
#cd /home/sskg8/code/DynamicShortestPath/PaRMAT/Release && make && ./PaRMAT -nVertices $vertices -nEdges $edges -sorted -noEdgeToSelf -noDuplicateEdges -threads 4 -output noWeightGraph.mtx && g++ -std=c++17 noWeightedtoWeightedMtx.cpp  -o noWeightedtoWeightedMtx && ./noWeightedtoWeightedMtx -i noWeightGraph.mtx -o $inputGraph -s 1 -e 10 && cp $inputGraph /home/sskg8/code/DynamicShortestPath/tools && cd /home/sskg8/code/DynamicShortestPath/tools

# #generate random edges
g++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx

g++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges2.mtx


# # Command 3 to run sosp code. g for input graph, c for changed edges, s for source. 
#g++ -std=c++17 SOSP.cpp -o program && ./program -g $inputGraph -c changeEdges.mtx -s $source

# #run parallal sosp
# g++ -std=c++17 -fopenmp parallelSOSPwithLoopSelectionWithoutCritical.cpp -o parallelSOSPwithLoopSelectionWithoutCritical -g -lstdc++ -lgomp 2> /dev/null && ./parallelSOSPwithLoopSelectionWithoutCritical -g $inputGraph -c changeEdges.mtx -s $source -t $threads

# # To generate MOSP
#g++ -std=c++17 -fopenmp MOSP.cpp -o MOSP -g -lstdc++ -lgomp 2> /dev/null  && ./MOSP -g1 $inputGraph  -c1 changeEdges.mtx -g2 $inputGraph -c2 changeEdges2.mtx -s 1

# # To generate Parallel MOSP
g++ -std=c++17 -fopenmp MOSPParallel.cpp -o MOSPParallel -g -lstdc++ -lgomp 2> /dev/null  && ./MOSPParallel -g1 $inputGraph  -c1 changeEdges.mtx -g2 $inputGraph -c2 changeEdges2.mtx -s 1 -t 4
 
## Vary threads
# for ((thread=1; thread<=$threads; thread+=4))
# do
#     g++ -std=c++17 -fopenmp parallelSOSPwithLoopSelectionWithoutCritical.cpp -o parallelSOSPwithLoopSelectionWithoutCritical -g -lstdc++ -lgomp 2> /dev/null && ./parallelSOSPwithLoopSelectionWithoutCritical -g $inputGraph -c changeEdges.mtx -s $source -t $thread >> performance.dat
# done

## Vary change edges
# rm "$dataset""_performance.dat"
# g++ -std=c++17  mtxDatasetToRandomWeight.cpp -o mtxDatasetToRandomWeight && ./mtxDatasetToRandomWeight -i $dataset -o $inputGraph -s 1 -e 10
# for ((ce=50; ce<=50000; ce*=5))
# do
#     echo -n $ce >> "$dataset""_performance.dat"
#     source='1'
#     insertedEdges=$ce
#     deletedEdges=$ce
#     inputGraph="weightedGraph.mtx"
#     dataset="inf-roadNet-CA.mtx"
#     threads="4"
#     g++ -std=c++17  mtxDatasetToRandomWeight.cpp -o mtxDatasetToRandomWeight && ./mtxDatasetToRandomWeight -i $dataset -o $inputGraph -s 1 -e 10
#     g++ -std=c++17 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i $insertedEdges -d $deletedEdges -g $inputGraph -o changeEdges.mtx
#     g++ -std=c++17 -fopenmp parallelSOSPwithLoopSelectionWithoutCritical.cpp -o parallelSOSPwithLoopSelectionWithoutCritical -g -lstdc++ -lgomp 2> /dev/null && ./parallelSOSPwithLoopSelectionWithoutCritical -g $inputGraph -c changeEdges.mtx -s $source -t $threads >> "$dataset""_performance.dat"
# done
