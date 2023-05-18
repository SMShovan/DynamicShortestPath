#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <limits>
#include "csr_conversion.h"

// CUDA kernel for Dijkstra's algorithm
__global__ void dijkstraKernel(const int* rowPointers, const int* columnIndices, const double* values,
                               double* distances, int* predecessors, bool* visited, int source, int n) {
    int u = blockIdx.x * blockDim.x + threadIdx.x;

    if (u < n) {
        distances[u] = (u == source) ? 0.0 : std::numeric_limits<double>::infinity();
        visited[u] = false;

        while (true) {
            int minVertex = -1;
            double minDistance = std::numeric_limits<double>::infinity();

            // Find the vertex with the minimum distance
            for (int v = 0; v < n; v++) {
                if (!visited[v] && distances[v] < minDistance) {
                    minVertex = v;
                    minDistance = distances[v];
                }
            }

            if (minVertex == -1) {
                break; // All vertices have been visited
            }

            visited[minVertex] = true;

            // Relax edges
            for (int i = rowPointers[minVertex]; i < rowPointers[minVertex + 1]; i++) {
                int v = columnIndices[i];
                double weight = values[i];

                if (!visited[v] && distances[minVertex] + weight < distances[v]) {
                    distances[v] = distances[minVertex] + weight;
                    predecessors[v] = minVertex;
                }
            }
        }
    }
}



std::vector<std::vector<Edge>> dijkstraParallelCUDA(const CSRMatrix& graph, int source) {
    int n = graph.rowPointers.size() - 1;

    // Allocate memory on the host
    double* h_distances = new double[n];
    int* h_predecessors = new int[n];
    bool* h_visited = new bool[n];

    // Allocate memory on the device
    double* d_distances;
    int* d_predecessors;
    bool* d_visited;
    cudaMalloc((void**)&d_distances, n * sizeof(double));
    cudaMalloc((void**)&d_predecessors, n * sizeof(int));
    cudaMalloc((void**)&d_visited, n * sizeof(bool));

    // Copy data from host to device
    cudaMemcpy(d_distances, h_distances, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_predecessors, h_predecessors, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_visited, h_visited, n * sizeof(bool), cudaMemcpyHostToDevice);

    // Launch the CUDA kernel
    int threadsPerBlock = 256;
    int numBlocks = (n + threadsPerBlock - 1) / threadsPerBlock;
    dijkstraKernel<<<numBlocks, threadsPerBlock>>>(graph.rowPointers.data(), graph.columnIndices.data(),
                                                   graph.values.data(), d_distances, d_predecessors,
                                                   d_visited, source, n);

    // Copy data from device to host
    cudaMemcpy(h_distances, d_distances, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_predecessors, d_predecessors, n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_visited, d_visited, n * sizeof(bool), cudaMemcpyDeviceToHost);

    // Free memory on the device
    cudaFree(d_distances);
    cudaFree(d_predecessors);
    cudaFree(d_visited);

    // Construct the single-source shortest path tree
    std::vector<std::vector<Edge>> shortestPaths(n);
    for (int i = 0; i < n; i++) {
        if (h_distances[i] != std::numeric_limits<double>::infinity()) {
            int currentNode = i;
            while (currentNode != source) {
                int predecessor = h_predecessors[currentNode];
                double weight = 0.0; // Default weight when edge weight is not available
                for (int j = graph.rowPointers[predecessor]; j < graph.rowPointers[predecessor + 1]; j++) {
                    if (graph.columnIndices[j] == currentNode) {
                        weight = graph.values[j];
                        break;
                    }
                }
                shortestPaths[predecessor].push_back(Edge(currentNode, weight));
                currentNode = predecessor;
            }
        }
    }

    // Free memory on the host
    delete[] h_distances;
    delete[] h_predecessors;
    delete[] h_visited;

    return shortestPaths;
}

int main(int argc, char** argv) {
    std::string inputFilename;
    int sourceNode = -1;

    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            inputFilename = argv[i + 1];
            i++;
        } else if (arg == "-s" && i + 1 < argc) {
            sourceNode = std::stoi(argv[i + 1]);
            i++;
        } else {
            std::cerr << "Invalid command-line argument: " << arg << std::endl;
            return 1;
        }
    }

    // Check if both -i and -s arguments are provided
    if (inputFilename.empty() || sourceNode == -1) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file> -s <source_node>" << std::endl;
        return 1;
    }

    // Convert the Matrix Market file to CSR format
    CSRMatrix csrMatrix = convertToCSR(inputFilename);

    // Perform Dijkstra's algorithm in parallel using CUDA
    std::vector<std::vector<Edge>> shortestPaths = dijkstraParallelCUDA(csrMatrix, sourceNode);

    // Print the CSR matrix representation of the shortest path tree
    printCSR(shortestPaths);

    return 0;
}