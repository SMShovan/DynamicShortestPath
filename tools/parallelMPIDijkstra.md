#Parallel Dijkstra's Algorithm with MPI

This program reads a graph from a Matrix Market file in CSR format, performs Dijkstra's algorithm for single-source shortest path in parallel using MPI, and prints the shortest distances from the source node to all other nodes.

#Requirements

* C++ compiler with MPI support
* Matrix Market file in CSR format

#Usage

Run the program with the following command-line arguments:

`mpirun -np <num_processes> ./dijkstraParallelMPI -i <input_file> -s <source_node>`

* `-i or --input`: Specifies the input file in Matrix Market format containing the graph data in CSR format.
* `-s or --source`: Specifies the source node from which to find the shortest paths.

#Example

Consider the following Matrix Market file graph.mtx representing a weighted graph in CSR format:

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

To find the single-source shortest paths from source node 1 in parallel using MPI, run the program as follows:

`mpirun -np 4 ./dijkstraParallelMPI -i graph.mtx -s 1`

The program will output the shortest distances from the source node to all other nodes:

```plaintext
Shortest Paths from Node 1:
Node 0: 0
Node 1: 0
Node 2: 3.7
Node 3: 4.2
Node 4: 6.1
```

The output displays the shortest distance from the source node (Node 1) to each node in the graph.

#Implementation Details

The program is implemented in C++ and utilizes MPI for parallelization. The main components of the program include:

* `csr_conversion.h`: A header file containing functions to convert a Matrix Market file to CSR format.
* `dijkstraParallelMPI.cpp`: The main C++ source code file that implements Dijkstra's algorithm for single-source shortest path in parallel using MPI.

The program follows these steps:

1. Each process initializes MPI and retrieves its rank and size.
1. Process 0 parses the command-line arguments to retrieve the input filename and the source node.
1. Process 0 converts the Matrix Market file to CSR format using the convertToCSR function.
1. Process 0 broadcasts the source node to all other processes.
1. Process 0 broadcasts the dimensions of the CSR matrix (number of rows, columns, and non-zeros) to all other processes.
1. All processes allocate memory for the CSR matrix based on the received dimensions.
1. Process 0 scatters the CSR matrix data (row pointers, column indices, and values) to all other processes using the MPI_Scatter function.
1. All processes create the adjacency list representation of the graph using the received CSR matrix data.
1. All processes call the dijkstraParallelMPI function to perform Dijkstra's algorithm in parallel. The function returns the shortest distances from the source node to all other nodes.
1. Process 0 gathers the shortest distances from all processes to construct the allDistances vector using the MPI_Gather function.
2. On the root process (rank 0), the printShortestPaths function is called to print the shortest distances for each node.
The program finalizes MPI using MPI_Finalize and returns 0.

These steps outline the main workflow of the program, including the MPI communication operations to distribute the data, perform the parallel Dijkstra's algorithm, and gather the results for printing on the root process.

Please note that proper error handling and additional MPI calls may be necessary to handle exceptional cases and ensure correct execution in real-world scenarios.