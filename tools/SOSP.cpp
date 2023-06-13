#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>

struct Edge {
    int source;
    int destination;
    double weight;

    Edge(int src, int dest, double w) : source(src), destination(dest), weight(w) {}
};

void printShortestPathTree(const std::vector<std::vector<int>>& ssspTree) {
    std::cout << "Shortest Path Tree:\n";

    for (int i = 0; i < ssspTree.size(); i++) {
        std::cout << "Node " << i + 1 << ": ";

        if (ssspTree[i].empty()) {
            std::cout << "No child nodes\n";
        } else {
            for (int child : ssspTree[i]) {
                std::cout << child << " ";
            }
            std::cout << "\n";
        }
    }
}


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

std::vector<std::vector<int>> dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode) {
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

    // Build the shortest path tree based on the shortest distances
    std::vector<std::vector<int>> ssspTree(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] != std::numeric_limits<double>::max()) {
            int parent = i + 1;
            for (const Edge& edge : graphCSR[i]) {
                int child = edge.destination + 1;
                if (shortestDist[child - 1] == shortestDist[i] + edge.weight) {
                    ssspTree[parent - 1].push_back(child);
                }
            }
        }
    }

    return ssspTree;
}

void updateShortestPath(std::vector<std::vector<int>>& ssspTree, const std::vector<std::vector<Edge>>& changedEdgesCSR) {
    for (const std::vector<Edge>& row : changedEdgesCSR) {
        for (const Edge& edge : row) {
            int sourceNode = edge.source;
            int destinationNode = edge.destination;
            double edgeWeight = edge.weight;

            // Check if both source and destination nodes are affected by the modification
            if (!ssspTree[sourceNode - 1].empty() && !ssspTree[destinationNode - 1].empty()) {
                // Update the shortest path tree for the affected vertices
                ssspTree[sourceNode - 1].push_back(destinationNode);
                ssspTree[destinationNode - 1].push_back(sourceNode);
            }
            // Check if only the source node is affected by the modification
            else if (!ssspTree[sourceNode - 1].empty()) {
                // Update the shortest path tree for the affected vertex
                ssspTree[sourceNode - 1].push_back(destinationNode);
            }
            // Check if only the destination node is affected by the modification
            else if (!ssspTree[destinationNode - 1].empty()) {
                // Update the shortest path tree for the affected vertex
                ssspTree[destinationNode - 1].push_back(sourceNode);
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

    // Find the initial shortest path tree using Dijkstra's algorithm
    std::vector<std::vector<int>> ssspTree = dijkstra(graphCSR, sourceNode);

    printShortestPathTree(ssspTree);

    // Read the changed edges file in Matrix Market format
    std::ifstream changedEdgesInputFile(changedEdgesFile);
    if (!changedEdgesInputFile.is_open()) {
        std::cerr << "Unable to open changed edges file: " << changedEdgesFile << std::endl;
        return 1;
    }

    // Convert the changed edges to CSR format
    std::vector<std::vector<Edge>> changedEdgesCSR = convertToCSR(changedEdgesInputFile);
    changedEdgesInputFile.close();

    // Update the shortest path tree based on the changed edges
    updateShortestPath(ssspTree, changedEdgesCSR);

    
    printShortestPathTree(ssspTree);
    

    return 0;
}
