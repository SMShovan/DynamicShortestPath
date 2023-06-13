#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

struct Edge {
    int source;
    int destination;
    double weight;

    Edge(int src, int dest, double w) : source(src), destination(dest), weight(w) {}
};

std::vector<std::vector<Edge>> convertToCSR(std::ifstream& inputFile) {
    std::string line;
    int numRows, numCols, numNonZero;
    if (!std::getline(inputFile, line)) {
        std::cerr << "Error: Unable to read the header line from the input file." << std::endl;
        exit(1);
    }
    std::istringstream headerStream(line);
    if (!(headerStream >> numRows >> numCols >> numNonZero)) {
        std::cerr << "Error: Invalid header format in the input file." << std::endl;
        exit(1);
    }

    std::vector<std::vector<Edge>> csrMatrix(numRows);

    int lineCount = 1;
    while (std::getline(inputFile, line)) {
        lineCount++;
        std::istringstream lineStream(line);
        int row, col;
        double value;
        if (!(lineStream >> row >> col >> value)) {
            std::cerr << "Error: Invalid line format at line " << lineCount << " in the input file." << std::endl;
            exit(1);
        }
        if (row < 1 || row > numRows || col < 1 || col > numCols) {
            std::cerr << "Error: Invalid vertex indices at line " << lineCount << " in the input file." << std::endl;
            exit(1);
        }
        csrMatrix[row - 1].emplace_back(row - 1, col - 1, value);
    }

    return csrMatrix;
}


std::vector<double> dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode) {
    int numNodes = graphCSR.size();

    // Create a vector to store the shortest distance from the source node to each node
    std::vector<double> shortestDist(numNodes, std::numeric_limits<double>::max());

    // Create a vector to track if a node has been visited during the algorithm
    std::vector<bool> visited(numNodes, false);

    // Set the distance of the source node to itself as 0
    shortestDist[sourceNode - 1] = 0;

    // Dijkstra's algorithm
    for (int i = 0; i < numNodes - 1; ++i) {
        // Find the node with the minimum distance among the unvisited nodes
        int minDistNode = -1;
        double minDist = std::numeric_limits<double>::max();
        for (int j = 0; j < numNodes; ++j) {
            if (!visited[j] && shortestDist[j] < minDist) {
                minDist = shortestDist[j];
                minDistNode = j;
            }
        }

        // Mark the minimum distance node as visited
        visited[minDistNode] = true;

        // Update the distances of the neighboring nodes
        for (const Edge& edge : graphCSR[minDistNode]) {
            int neighbor = edge.destination;
            double weight = edge.weight;
            if (!visited[neighbor] && shortestDist[minDistNode] + weight < shortestDist[neighbor]) {
                shortestDist[neighbor] = shortestDist[minDistNode] + weight;
            }
        }
    }

    return shortestDist;
}


void updateShortestPath(std::vector<double>& shortestPath, const std::vector<std::vector<Edge>>& changedEdgesCSR) {
    // Iterate over the rows in the changed edges CSR matrix
    for (const std::vector<Edge>& row : changedEdgesCSR) {
        for (const Edge& edge : row) {
            int sourceNode = edge.source;
            int destinationNode = edge.destination;
            double edgeWeight = edge.weight;

            // Check if both source and destination nodes are affected by the modification
            if (shortestPath[sourceNode] != std::numeric_limits<double>::max() &&
                shortestPath[destinationNode] != std::numeric_limits<double>::max()) {
                // Update the CSR matrix for the affected vertices
                shortestPath[sourceNode] = edgeWeight;
                shortestPath[destinationNode] = edgeWeight;
            }
            // Check if only the source node is affected by the modification
            else if (shortestPath[sourceNode] != std::numeric_limits<double>::max()) {
                // Update the CSR matrix for the affected vertex
                shortestPath[sourceNode] = edgeWeight;
            }
            // Check if only the destination node is affected by the modification
            else if (shortestPath[destinationNode] != std::numeric_limits<double>::max()) {
                // Update the CSR matrix for the affected vertex
                shortestPath[destinationNode] = edgeWeight;
            }
            // Neither source nor destination node is affected, no changes are required
        }
    }
}


int main(int argc, char** argv) {
    // Check the command-line arguments
    if (argc < 4) {
        std::cerr << "Usage: ./program -g <graph_file> -c <changed_edges_file> -s <source_node>\n";
        return 1;
    }

    std::string graphFile;
    std::string changedEdgesFile;
    int sourceNode;

    // Parse the command-line arguments
    for (int i = 1; i < argc; i += 2) {
        std::string option(argv[i]);
        std::string argument(argv[i + 1]);

        if (option == "-g") {
            graphFile = argument;
        } else if (option == "-c") {
            changedEdgesFile = argument;
        } else if (option == "-s") {
            sourceNode = std::stoi(argument);
        } else {
            std::cerr << "Invalid option: " << option << "\n";
            return 1;
        }
    }

    std::ifstream graphInputFile(graphFile);
    if (!graphInputFile.is_open()) {
        std::cerr << "Unable to open graph file: " << graphFile << std::endl;
        return 1;
    }

    // Convert the input graph to CSR format
    std::vector<std::vector<Edge>> graphCSR = convertToCSR(graphInputFile);
    graphInputFile.close();

    // Find the initial shortest path using Dijkstra's algorithm
    std::vector<double> shortestPath = dijkstra(graphCSR, sourceNode);

    // Read the changed edges file in Matrix Market format
    std::ifstream changedEdgesInputFile(changedEdgesFile);
    if (!changedEdgesInputFile.is_open()) {
        std::cerr << "Unable to open changed edges file: " << changedEdgesFile << std::endl;
        return 1;
    }

    // Convert the changed edges to CSR format
    std::vector<std::vector<Edge>> changedEdgesCSR = convertToCSR(changedEdgesInputFile);
    changedEdgesInputFile.close();

    // Update the shortest path based on the changed edges
    updateShortestPath(shortestPath, changedEdgesCSR);

    // Print or store the updated shortest path
    std::cout << "Updated Shortest Path from Node " << sourceNode << ":" << std::endl;
    for (int i = 0; i < shortestPath.size(); i++) {
        std::cout << "Node " << i + 1 << ": " << shortestPath[i] << std::endl;
    }

    return 0;
}


