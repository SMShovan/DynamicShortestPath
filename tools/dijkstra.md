# Dijkstra's Algorithm for Single-Source Shortest Path

This program reads a graph from a Matrix Market file in CSR format, performs Dijkstra's algorithm for single-source shortest path, and prints the CSR matrix representation of the resulting shortest path tree.

## Requirements

- C++ compiler
- Matrix Market file in CSR format

## Usage

Run the program with the following command-line arguments:

`./dijkstra -i <input_file> -s <source_node>`


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
To find the single-source shortest path tree from source node 1, run the program as follows:

`./dijkstra -i graph.mtx -s 1`

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
* `dijkstra.cpp`: The main C++ source code file that implements Dijkstra's algorithm for single-source shortest path.

The program follows these steps:

1. Parse the command-line arguments to retrieve the input filename and the source node.
1. Convert the Matrix Market file to CSR format using the convertToCSR function.
1. Perform Dijkstra's algorithm using the dijkstra function, which returns the single-source shortest path tree as a CSR matrix.
1. Print the CSR matrix representation of the shortest path tree using the printCSR function.

#Limitations

The program assumes that the input graph is represented in Matrix Market format and follows the CSR format specification.
The input graph should be a weighted graph. If the Matrix Market file contains only the pattern of the graph, the program assumes that all edges have a weight of 1.0.
Dependencies

* C++ Standard Library

