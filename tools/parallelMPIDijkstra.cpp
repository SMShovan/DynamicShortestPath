#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <limits>
#include "csr_conversion.h"
#include <mpi.h>

struct Edge {
    int destination;
    double weight;

    Edge(int dest, double w) : destination(dest), weight(w) {}
};

std::vector<double> dijkstraParallelMPI(const std::vector<std::vector<Edge>>& graph, int sourceNode) {
    int numNodes = graph.size();

    // Initialize distance array with infinity
    std::vector<double> distance(numNodes, std::numeric_limits<double>::infinity());
    distance[sourceNode] = 0;

    // Create a priority queue to store nodes
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    pq.push(std::make_pair(0, sourceNode));

    // Process the nodes until the priority queue becomes empty
    while (!pq.empty()) {
        int currentNode = pq.top().second;
        double currentDistance = pq.top().first;
        pq.pop();

        // Skip if the node has already been processed
        if (currentDistance > distance[currentNode]) {
            continue;
        }

        // Traverse the adjacent nodes of the current node
        for (const auto& edge : graph[currentNode]) {
            int neighborNode = edge.destination;
            double weight = edge.weight;
            double newDistance = currentDistance + weight;

            // Relax the edge if a shorter path is found
            if (newDistance < distance[neighborNode]) {
                distance[neighborNode] = newDistance;
                pq.push(std::make_pair(newDistance, neighborNode));
            }
        }
    }

    return distance;
}

void printShortestPaths(const std::vector<double>& distances, int sourceNode) {
    int numNodes = distances.size();

    std::cout << "Shortest Paths from Node " << sourceNode << ":\n";
    for (int i = 0; i < numNodes; ++i) {
        std::cout << "Node " << i << ": ";
        if (distances[i] == std::numeric_limits<double>::infinity()) {
            std::cout << "Unreachable";
        } else {
            std::cout << distances[i];
        }
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string inputFilename;
    int sourceNode = -1;

    if (rank == 0) {
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
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Check if both -i and -s arguments are provided
        if (inputFilename.empty() || sourceNode == -1) {
            std::cerr << "Usage: " << argv[0] << " -i <input_file> -s <source_node>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast the source node to all processes
    MPI_Bcast(&sourceNode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Convert the Matrix Market file to CSR format on the root process
    CSRMatrix csrMatrix;
    if (rank == 0) {
        csrMatrix = convertToCSR(inputFilename);
    }

    // Broadcast the CSR matrix dimensions to all processes
    int numRows, numCols, numNonZeros;
    if (rank == 0) {
        numRows = csrMatrix.rowPointers.size() - 1;
        numCols = csrMatrix.columnIndices.size();
        numNonZeros = csrMatrix.values.size();
    }
    MPI_Bcast(&numRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numNonZeros, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for the CSR matrix on each process
    std::vector<int> rowPointers(numRows + 1);
    std::vector<int> columnIndices(numNonZeros);
    std::vector<double> values(numNonZeros);

    // Scatter the CSR matrix data to all processes
    MPI_Scatter(csrMatrix.rowPointers.data(), numRows + 1, MPI_INT, rowPointers.data(), numRows + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(csrMatrix.columnIndices.data(), numNonZeros, MPI_INT, columnIndices.data(), numNonZeros, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(csrMatrix.values.data(), numNonZeros, MPI_DOUBLE, values.data(), numNonZeros, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Create the adjacency list representation of the graph
    std::vector<std::vector<Edge>> graph(numRows);
    for (int i = 0; i < numRows; ++i) {
        for (int j = rowPointers[i]; j < rowPointers[i + 1]; ++j) {
            int destination = columnIndices[j];
            double weight = values[j];
            graph[i].push_back(Edge(destination, weight));
        }
    }

    // Perform Dijkstra's algorithm in parallel using MPI
    std::vector<double> shortestDistances = dijkstraParallelMPI(graph, sourceNode);

    // Gather the shortest distances from all processes to the root process
    std::vector<double> allDistances(numRows);
    MPI_Gather(shortestDistances.data(), numRows, MPI_DOUBLE, allDistances.data(), numRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print the shortest distances on the root process
    if (rank == 0) {
        printShortestPaths(allDistances, sourceNode);
    }

    MPI_Finalize();

    return 0;
}
