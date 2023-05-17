# Parallel Dijkstra's Algorithm for Single-Source Shortest Path

This program reads a graph from a Matrix Market file in CSR format, performs Dijkstra's algorithm for single-source shortest path in parallel using OpenMP, and prints the CSR matrix representation of the resulting shortest path tree.

## Requirements

- C++ compiler with OpenMP support
- Matrix Market file in CSR format

## Usage

Run the program with the following command-line arguments:

`./parallelDijkstra -i <input_file> -s <source_node>`


- `-i` or `--input`: Specifies the input file in Matrix Market format containing the graph data in CSR format.
- `-s` or `--source`: Specifies the source node from which to find the shortest paths.

## Example

Consider the following Matrix Market file `graph.mtx` representing a weighted graph in CSR format:

```plaintext
%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This is a comment line
%
%=================================================================================
5 5 8
1 1 2.0
1 3 1.5
2 2 3.7
3 1 4.2
3 4 2.8
4 3 6.1
5 2 5.0
5 5 9.3
```

To find the single-source shortest path tree from source node 1 in parallel, run the program as follows:

`./parallelDijkstra -i graph.mtx -s 1`

The program will output the CSR matrix representation of the resulting shortest path tree:

```plaintext
CSR Matrix Representation:
Row Pointers: 0 2 3 5 6 8
Column Indices: 0 2 0 3 1 3 2 4
Values: 2 1.5 3.7 4.2 2.8 6.1 5 9.3
```

The row pointers indicate the starting indices of each row in the column indices and values arrays. The column indices represent the adjacent nodes of each row, and the values represent the corresponding edge weights.

#Implementation Details

The program is implemented in C++ and consists of the following main components:

* `csr_conversion.h`: A header file containing functions to convert a Matrix Market file to CSR format.

* `parallelDijkstra.cpp`: The main C++ source code file that implements Dijkstra's algorithm for single-source shortest path in parallel using OpenMP.

The program follows these steps:

1. Parse the command-line arguments to retrieve the input filename and the source node.
1. Convert the Matrix Market file to CSR format using the convertToCSR function.
1. Perform Dijkstra's algorithm in parallel using OpenMP by calling the dijkstraParallel function, which returns the single-source shortest path tree as a CSR matrix.
1. Print the CSR matrix representation of the shortest path tree using the printCSR function.

#Parallelization

The Dijkstra's algorithm is parallelized using OpenMP to take advantage of multiple threads and expedite the computation. The dijkstraParallel function parallelizes the traversal of adjacent nodes of each node using a #pragma omp parallel for directive, ensuring efficient exploration of the CSR data structure.

# Utilizing CSR Format for Efficiency

The CSR format is a popular choice for representing sparse matrices due to its memory efficiency and efficient traversal of adjacent nodes. In the parallelDijkstra.cpp code, the CSR format is utilized to improve the efficiency of the Dijkstra's algorithm.

1. Construction of CSR Matrix: The program starts by converting the input Matrix Market file to the CSR format using the convertToCSR function. The CSR format stores the graph's non-zero values and their corresponding row and column indices in three separate arrays: rowPointers, columnIndices, and values. This format allows for efficient access to the graph's structure and weights.
1. Traversal of Adjacent Nodes: In the dijkstraParallel function, the parallelization of the Dijkstra's algorithm takes advantage of the CSR format. The traversal of adjacent nodes for each node is performed using a #pragma omp parallel for directive. This directive parallelizes the loop, allowing multiple threads to explore the adjacent nodes concurrently.
	* 	Row Pointers: The rowPointers array of the CSR format is used to determine the starting index and ending index of the adjacent nodes for each node. This information enables efficient access to the adjacent nodes during parallel traversal.
	* 	Column Indices and Values: The columnIndices array provides the column indices of the adjacent nodes, and the values array contains the corresponding edge weights. These arrays are accessed to determine the destination nodes and weights during the relaxation process of Dijkstra's algorithm.

By utilizing the CSR format, the parallelDijkstra.cpp code avoids unnecessary memory overhead and ensures efficient traversal of adjacent nodes during the parallelized Dijkstra's algorithm. The CSR format's structure enables efficient access to the graph data, resulting in improved performance and reduced computational complexity for large sparse graphs.

This efficient utilization of the CSR format helps the parallelDijkstra.cpp program handle large graphs more effectively, making it suitable for applications where efficiency and scalability are crucial factors.


#Limitations

The program assumes that the input graph is represented in Matrix Market format and follows the CSR format specification.
The input graph should be a weighted graph. If the Matrix Market file contains only the pattern of the graph, the program assumes that all edges have a weight of 1.0.

#Dependencies

C++