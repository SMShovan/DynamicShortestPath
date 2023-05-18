# Parallel Dijkstra's Algorithm for Single-Source Shortest Path (CUDA)

This program reads a graph from a Matrix Market file in CSR format, performs Dijkstra's algorithm for single-source shortest path in parallel using CUDA GPU, and prints the CSR matrix representation of the resulting shortest path tree.

## Requirements

- C++ compiler with CUDA support
- CUDA-enabled GPU
- Matrix Market file in CSR format

## Usage

Run the program with the following command-line arguments:

`./dijkstraParallelCUDA -i <input_file> -s <source_node>`


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

To find the single-source shortest path tree from source node 1 in parallel using CUDA, run the program as follows:

`./dijkstraParallelCUDA -i graph.mtx -s 1`

The program will output the CSR matrix representation of the resulting shortest path tree:

```plaintext
CSR Matrix Representation:
Row Pointers: 0 2 3 5 6 8
Column Indices: 0 2 0 3 1 3 2 4
Values: 2 1.5 3.7 4.2 2.8 6.1 5 9.3
```

The row pointers indicate the starting indices of each row in the column indices and values arrays. The column indices represent the adjacent nodes of each row, and the values represent the corresponding edge weights.

#Implementation Details

The program is implemented in C++ and utilizes CUDA GPU for parallelization. The main components of the program include:

* `csr_conversion.h`: A header file containing functions to convert a Matrix Market file to CSR format.
* `dijkstraParallelCUDA.cpp`: The main C++ source code file that implements Dijkstra's algorithm for single-source shortest path in parallel using CUDA.

The program follows these steps:

1. Parse the command-line arguments to retrieve the input filename and the source node.
1. Convert the Matrix Market file to CSR format using the convertToCSR function.
1. Perform Dijkstra's algorithm in parallel using CUDA GPU by calling the dijkstraParallelCUDA function, which returns the single-source shortest path tree as a CSR matrix.
1. Print the CSR matrix representation of the shortest path tree using the printCSR function.

#Utilizing CUDA GPU for Parallelization

The program utilizes CUDA GPU for parallelization, enabling faster computation of Dijkstra's algorithm. The key steps in leveraging CUDA GPU for parallel execution are as follows:

1. CUDA Kernel: The dijkstraKernel is a CUDA kernel function that performs Dijkstra's algorithm in parallel on the GPU. It utilizes CUDA parallelism by launching multiple threads to process different nodes concurrently. Each thread processes a subset of nodes using a combination of thread indexing and block indexing.
1. Memory Management: The program allocates memory on both the host and device. The input graph data in CSR format is stored on the host, while the intermediate results and output data are stored on the device. The cudaMalloc function is used to allocate device memory, and the cudaMemcpy function is used to transfer data between the host and device.
1. Parallel Processing: The CUDA kernel is launched with a specified number of blocks and threads per block. Each thread is responsible for processing a subset of nodes in parallel, performing the necessary operations such as relaxation of edges and updating the shortest distance and predecessor information. The CUDA kernel exploits the parallel architecture of the GPU, allowing for efficient computation of the algorithm.
1. Data Transfer: After the CUDA kernel execution, the results are copied back from the device to the host using the cudaMemcpy function. The host memory is then used to construct the single-source shortest path tree representation.
1. Memory Cleanup: Memory allocated on the device is freed using the cudaFree function to avoid memory leaks.

By utilizing CUDA GPU for parallelization, the program achieves significant speedup compared to sequential execution. The parallel execution allows for simultaneous processing of multiple nodes, taking advantage of the GPU's parallel architecture and high computational power. This approach is especially beneficial for large-scale graphs, where the parallelism provided by the GPU can greatly reduce the overall computation time.
