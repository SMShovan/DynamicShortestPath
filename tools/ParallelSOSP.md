# Parallel Shared Memory Shortest Path Tree (SOSP) Algorithm

## Introduction
The Parallel Shared Memory Shortest Path Tree (SOSP) algorithm is an optimized version of the SOSP algorithm that utilizes parallelism using OpenMP for efficient computation of the shortest path tree in a graph. This algorithm is suitable for shared memory systems with multiple processing units.

## Usage
The code is implemented in C++ and requires the OpenMP library. To compile and run the code, follow these steps:

1. Save the code in a file named `parallelSharedMemorySOSP.cpp`.
2. Compile the code using a C++ compiler with OpenMP support. For example:
```
g++ -fopenmp parallelSharedMemorySOSP.cpp -o parallelSharedMemorySOSP
```
3. Run the executable:

```
./parallelSharedMemorySOSP -g <graph_file> -c <changed_edges_file> -s <source_node>
```
Replace `<graph_file>` with the path to the graph file in Matrix Market format, `<changed_edges_file>` with the path to the file containing changed edges in Matrix Market format, and `<source_node>` with the index of the source node.

## Algorithm Description
The parallelSharedMemorySOSP algorithm consists of the following steps:

1. Input Parsing: The program reads the input graph and changed edges files in Matrix Market format and converts them to Compressed Sparse Row (CSR) format.
2. Dijkstra's Algorithm: The algorithm performs Dijkstra's algorithm to compute the initial shortest path tree from the source node. This step is parallelized using OpenMP to accelerate the computation.
3. Shortest Path Tree Update: The algorithm updates the shortest path tree based on the changed edges. It identifies inserted and deleted edges, modifies the shortest distances, and updates the parent-child relationships in the tree. This step is also parallelized using OpenMP for efficient processing.
4. Output: The final shortest path tree is printed to the console.

## Functions

### `convertToCSR`
This function converts the input graph in Matrix Market format to Compressed Sparse Row (CSR) format. It takes an input file stream as a parameter and returns a vector of vectors representing the CSR matrix.

```cpp
std::vector<std::vector<Edge>> convertToCSR(std::ifstream& inputFile);
```
###Dijkstra
This function performs Dijkstra's algorithm to compute the initial shortest path tree from the source node. It takes the graph in CSR format, the source node index, and a vector to store the shortest distances as parameters. It returns a vector of vectors representing the shortest path tree.
```
std::vector<std::vector<int>> dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode, std::vector<double>& shortestDist);
```
###markSubtreeAffected
This function marks the subtree rooted at a given node as affected. It takes the parent-child representation of the tree, the shortest distances, a vector to track affected nodes, a vector to track affected and deleted nodes, and the node index as parameters.

```
void markSubtreeAffected(std::vector<std::pair<int, std::vector<int>>>& parentChildSSP, std::vector<double>& shortestDist, std::vector<bool>& affectedNodes, std::vector<bool>& affectedDelNodes, int node);
```

###convertToParentChildSSP
This function converts the shortest path tree in vector of vectors representation to parent-child representation. It takes the vector of vectors representing the tree as a parameter and returns a vector of pairs representing the parent-child relationships.

```
std::vector<std::pair<int, std::vector<int>>> convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree);
```

###updateShortestPath
This function updates the shortest path tree based on the changed edges. It takes the parent-child representation of the tree, the graph in CSR format, the changed edges in CSR format, and the shortest distances as parameters.
```
void updateShortestPath(std::vector<std::pair<int, std::vector<int>>>& ssspTree, const std::vector<std::vector<Edge>>& graphCSR, const std::vector<std::vector<Edge>>& changedEdgesCSR, std::vector<double>& shortestDist);
```
##Complexity Analysis

The time complexity of the parallelSharedMemorySOSP algorithm depends on the size of the input graph and the number of iterations required for Dijkstra's algorithm. The overall complexity can be analyzed as follows:

* Converting the input graph to CSR format: O(E + V), where E is the number of edges and V is the number of vertices.
* Dijkstra's algorithm: O(V^2) in the worst case.
* Shortest path tree update: O(E + V) for identifying inserted and deleted edges, and O(E) for updating the shortest distances and parent-child relationships.
* Printing the shortest path tree: O(V).
* Therefore, the overall time complexity of the algorithm is O(V^2 + E), where V is the number of vertices and E is the number of edges in the input graph.

The space complexity of the algorithm is O(V^2) to store the graph in CSR format and the shortest path tree.

##Example

Here's an example of how to run the code:
```
./parallelSharedMemorySOSP -g graph.txt -c changed_edges.txt -s 1
```
In this example, the graph is stored in the file graph.txt, the changed edges are stored in the file changed_edges.txt, and the source node index is 1.

The program will compute the shortest path tree from the source node and print the resulting tree to the console.

```plaintext
Shortest Path Tree:
Node 1: 2 3 
Node 2: 4 5 
Node 4: 6 
Node 6:
Node 3:
Node 5:
```
This output indicates that the source node is connected to nodes 2 and 3, node 2 is connected to nodes 4 and 5, node 4 is connected to node 6, and nodes 3 and 5 have no child nodes.


## Conclusion
The parallelSharedMemorySOSP algorithm provides an efficient parallel implementation of the SOSP algorithm using OpenMP. By leveraging multiple processing units in shared memory systems, it improves the computation speed of finding the shortest path tree in a graph. The algorithm can be used in various applications that require efficient graph traversal and shortest path calculations.

Feel free to use and modify the code according to your requirements. Happy coding!




