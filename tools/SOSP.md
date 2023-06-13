# Shortest Path Tree Update Documentation

This documentation provides an overview of the Shortest Path Tree Update program, including its purpose, usage, function descriptions, examples, and complexity analysis.

## Purpose

The Shortest Path Tree Update program is designed to update the shortest path tree of a graph based on changes to the graph's edges. It takes as input a graph file in Compressed Sparse Row (CSR) format, a file containing the changed edges in Matrix Market format, and a source node for the shortest path tree. The program then updates the shortest path tree based on the changes and outputs the updated tree.

## Usage

To run the Shortest Path Tree Update program, use the following command-line format:

```
./program -g <graph_file> -c <changed_edges_file> -s <source_node>
```


- `<graph_file>`: The path to the graph file in Compressed Sparse Row (CSR) format.
- `<changed_edges_file>`: The path to the file containing the changed edges in Matrix Market format.
- `<source_node>`: The source node for the shortest path tree.

## Function Descriptions

The program consists of several functions to convert the input graph to CSR format, perform Dijkstra's algorithm to find the initial shortest path tree, update the shortest path tree based on changed edges, and print the updated tree.

### `convertToCSR(std::ifstream& inputFile)`

This function takes an input file stream of the graph in Matrix Market format and converts it to Compressed Sparse Row (CSR) format. It returns a 2D vector representing the CSR matrix.

### `dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode, std::vector<double>& shortestDist)`

This function performs Dijkstra's algorithm on the graph represented in CSR format to find the shortest path tree. It takes the CSR graph, the source node, and a reference to a vector to store the shortest distances. It returns a 2D vector representing the shortest path tree.

### `printShortestPathTree(const std::vector<std::pair<int, std::vector<int>>>& parentChildSSP)`

This function prints the updated shortest path tree in parent-child format. It takes a vector of pairs representing the parent-child relationships in the tree.

### `markSubtreeAffected(std::vector<std::pair<int, std::vector<int>>>& parentChildSSP, std::vector<double>& shortestDist, std::vector<bool>& affectedNodes,std::vector<bool>& affectedDelNodes, int node)`

This function marks the subtree rooted at a given node as affected due to edge deletions. It updates the affected nodes vector, the shortest distances vector, and the parent-child shortest path tree vector.

### `convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree)`

This function converts the shortest path tree represented as a 2D vector to a parent-child relationship data structure. It takes the shortest path tree and returns a vector of pairs representing the parent-child relationships.

### `updateShortestPath(std::vector<std::pair<int, std::vector<int>>>& ssspTree, const std::vector<std::vector<Edge>>& graphCSR, const std::vector<std::vector<Edge>>& changedEdgesCSR, std::vector<double>& shortestDist)`

This function updates the shortest path tree based on the changed edges. It takes the parent-child shortest path tree, the original CSR graph, the changed edges CSR graph, and the shortest distances vector. It modifies the parent-child shortest path tree and shortest distances based on the changes.

## Examples

Here are some example usages of the Shortest Path Tree Update program:

1. Update the shortest path tree for a graph with changes:
```
./program -g graph.txt -c changes.txt -s 1
```

In this example, the program reads the graph from `graph.txt`, the changes from `changes.txt`, and sets node 1 as the source node for the shortest path tree.

2. Update the shortest path tree for a large graph:
```
./program -g graph_large.txt -c changes_large.txt -s 10
```


This example demonstrates the usage of the program with a large graph and a specific source node.

Feel free to modify the program and experiment with different graph files and changed edge files to update the shortest path tree according to your requirements.

## Complexity Analysis

The time complexity of the Shortest Path Tree Update program depends on the size of the input graph and the number of changed edges. The key operations are:

- Converting the input graph to CSR format: O(|E|), where |E| is the number of edges in the graph.
- Performing Dijkstra's algorithm: O((|V| + |E|) log |V|), where |V| is the number of vertices and |E| is the number of edges in the graph.
- Updating the shortest path tree based on changed edges: O(|E'| + |V|), where |E'| is the number of changed edges and |V| is the number of vertices in the graph.
- Printing the updated shortest path tree: O(|V| + |E|), where |V| is the number of vertices and |E| is the number of edges in the graph.

The space complexity of the program is O(|V| + |E|), where |V| is the number of vertices and |E| is the number of edges in the graph. This accounts for the storage of the CSR graph, changed edges CSR graph, shortest distances, parent-child shortest path tree, and other auxiliary data structures.

---

Feel free to copy the entire documentation, including the code, function descriptions, examples, and complexity analysis. If you have any further questions, please let me know!
