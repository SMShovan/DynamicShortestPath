#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <limits>
#include "csr_conversion.h"

// Structure to represent a graph edge
struct Edge {
    int destination;
    double weight;

    Edge(int dest, double w) : destination(dest), weight(w) {}
};

// Structure to represent a graph node
struct Node {
    int id;
    double distance;

    Node(int _id, double _dist) : id(_id), distance(_dist) {}

    // Overload the comparison operator for priority queue
    bool operator<(const Node& other) const {
        return distance > other.distance;
    }
};

// Function to perform Dijkstra's algorithm for single-source shortest path
std::vector<std::vector<Edge>> dijkstra(const CSRMatrix& graph, int source) {
    int n = graph.rowPointers.size() - 1;

    // Initialize distance array with infinity
    std::vector<double> distance(n, std::numeric_limits<double>::infinity());

    // Create a priority queue to store nodes
    std::priority_queue<Node> pq;

    // Set the distance of the source node to 0 and push it into the priority queue
    distance[source] = 0;
    pq.push(Node(source, 0));

    // Process the nodes until the priority queue becomes empty
    while (!pq.empty()) {
        // Get the node with the minimum distance from the priority queue
        Node current = pq.top();
        pq.pop();

        int u = current.id;
        double dist = current.distance;

        // Skip if the node has already been processed
        if (dist > distance[u]) {
            continue;
        }

        // Traverse the adjacent nodes of the current node
        for (int i = graph.rowPointers[u]; i < graph.rowPointers[u + 1]; i++) {
            int v = graph.columnIndices[i];
            double weight = graph.values[i];

            // Relax the edge if a shorter path is found
            if (dist + weight < distance[v]) {
                distance[v] = dist + weight;
                pq.push(Node(v, distance[v]));
            }
        }
    }

    // Construct and return the single-source shortest path tree as a CSR matrix
    std::vector<std::vector<Edge>> shortestPaths(n);
    for (int i = 0; i < n; i++) {
        for (int j = graph.rowPointers[i]; j < graph.rowPointers[i + 1]; j++) {
            int v = graph.columnIndices[j];
            double weight = graph.values[j];
            if (distance[v] == distance[i] + weight) {
                shortestPaths[i].push_back(Edge(v, weight));
            }
        }
    }

    return shortestPaths;
}

void printCSR(const std::vector<std::vector<Edge>>& csrMatrix) {
    int n = csrMatrix.size();

    std::cout << "CSR Matrix Representation:\n";
    std::cout << "Row Pointers: ";
    int rowPointer = 0;
    for (const auto& row : csrMatrix) {
        std::cout << rowPointer << " ";
        rowPointer += row.size();
    }
    std::cout << rowPointer << "\nColumn Indices: ";
    for (const auto& row : csrMatrix) {
        for (const auto& edge : row) {
            std::cout << edge.destination << " ";
        }
    }
    std::cout << "\nValues: ";
    for (const auto& row : csrMatrix) {
        for (const auto& edge : row) {
            std::cout << edge.weight << " ";
        }
    }
    std::cout << std::endl;
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

    // Perform Dijkstra's algorithm
    std::vector<std::vector<Edge>> shortestPaths = dijkstra(csrMatrix, sourceNode);

    // Print the CSR matrix representation of the single-source shortest path tree
    printCSR(shortestPaths);

    return 0;
}
