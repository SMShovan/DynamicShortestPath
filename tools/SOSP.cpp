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

std::vector<std::vector<int>> dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode, std::vector<double>& shortestDist) {
    int numNodes = graphCSR.size();

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



// void printShortestPathTree(const std::vector<std::vector<int>>& ssspTree) {
//     std::cout << "Shortest Path Tree:\n";

//     for (int i = 0; i < ssspTree.size(); i++) {
//         std::cout << "Node " << i + 1 << ": ";

//         if (ssspTree[i].empty()) {
//             std::cout << "No child nodes\n";
//         } else {
//             for (int child : ssspTree[i]) {
//                 std::cout << child << " ";
//             }
//             std::cout << "\n";
//         }
//     }
// }

void printShortestPathTree(const std::vector<std::pair<int, std::vector<int>>>& parentChildSSP) {
    std::cout << "Shortest Path Tree:\n";
    for (const auto& node : parentChildSSP) {
        std::cout << "Node " << node.first << ": ";
        for (int child : node.second) {
            std::cout << child << " ";
        }
        std::cout << "\n";
    }
}


void markSubtreeAffected(std::vector<std::pair<int, std::vector<int>>>& parentChildSSP, std::vector<double>& shortestDist, std::vector<bool>& affectedNodes,std::vector<bool>& affectedDelNodes, int node) {
    affectedNodes[node - 1] = true;
    affectedDelNodes[node - 1] = false;
    shortestDist[node - 1] = std::numeric_limits<double>::infinity();


    for (int child : parentChildSSP[node - 1].second) {
        markSubtreeAffected(parentChildSSP, shortestDist, affectedNodes, affectedDelNodes,  child);
    }
}


std::vector<std::pair<int, std::vector<int>>> convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree) {
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP(ssspTree.size());

    for (int i = 0; i < ssspTree.size(); ++i) {
        parentChildSSP[i].first = i + 1;  // Set the node as the parent in the parent-child SSP structure
        parentChildSSP[i].second = ssspTree[i];  // Set the child nodes
    }

    return parentChildSSP;
}



void updateShortestPath(std::vector<std::pair<int, std::vector<int>>>& ssspTree, const std::vector<std::vector<Edge>>& graphCSR, const std::vector<std::vector<Edge>>& changedEdgesCSR, std::vector<double>& shortestDist) {
    std::vector<Edge> insertedEdges;
    std::vector<Edge> deletedEdges;
    std::vector<bool> affectedNodes(ssspTree.size(), false);
    std::vector<bool> affectedDelNodes(ssspTree.size(), false);

    // Identify inserted and deleted edges
    for (const std::vector<Edge>& row : changedEdgesCSR) {
        for (const Edge& edge : row) {
            if (edge.weight > 0) {
                insertedEdges.push_back(edge);
            } else if (edge.weight < 0) {
                deletedEdges.push_back(edge);
            }
        }
    }

    // Print the inserted edges and identify affected nodes
    std::cout << "Inserted Edges:\n";
    for (const Edge& edge : insertedEdges) {
        std::cout << "Source: " << edge.source << " Destination: " << edge.destination << " Weight: " << edge.weight << "\n";
        
        int x,y; 
        if(shortestDist[edge.source - 1] > shortestDist[edge.destination - 1])
            x = edge.destination - 1; 
        else
            y = edge.source - 1; 
        
        if (shortestDist[y - 1] > shortestDist[x - 1]) {
            shortestDist[y - 1] = shortestDist[x - 1] + edge.weight;
            ssspTree[y - 1].first = x - 1; 
            affectedNodes[y - 1] = true; 
        }
    }

    // Print the deleted edges and mark affected nodes
    std::cout << "Deleted Edges:\n";
    for (const Edge& edge : deletedEdges) {
        std::cout << "Source: " << edge.source << " Destination: " << edge.destination << " Weight: " << edge.weight << "\n";
        affectedNodes[edge.destination - 1] = true;
        affectedDelNodes[edge.destination - 1] = true;
        shortestDist[edge.destination - 1] = std::numeric_limits<double>::infinity();
        markSubtreeAffected(ssspTree, shortestDist, affectedNodes, affectedDelNodes, edge.destination);
    }

    bool hasAffectedNodes = true;
    while (hasAffectedNodes) {
        hasAffectedNodes = false;

        // Check for affected nodes and update distances
        for (int v = 0; v < ssspTree.size(); ++v) {
            if (affectedNodes[v - 1]) {
                affectedNodes[v - 1] = false;  // Reset affected flag

                for (const Edge& edge : graphCSR[v - 1]) {
                    int n = edge.destination;

                    if (shortestDist[v - 1] > shortestDist[n - 1] + edge.weight) {
                        shortestDist[v - 1] = shortestDist[n - 1] + edge.weight;
                        ssspTree[v - 1].first = n;
                        affectedNodes[v - 1] = true;  // Mark as affected
                        hasAffectedNodes = true;  // Set flag for next iteration
                    }
                }
            }
        }
    }



    // Print the affected nodes
    // std::cout << "Affected Nodes:\n";
    // for (int i = 0; i < affectedNodes.size(); i++) {
    //     if (affectedNodes[i]) {
    //         std::cout << "Node " << i + 1 << "\n";
    //     }
    // }

    // Propagate

    


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
    std::vector<double> shortestDist(graphCSR.size(), std::numeric_limits<double>::max());
    std::vector<std::vector<int>> ssspTree = dijkstra(graphCSR, sourceNode, shortestDist);

    // Convert the ssspTree to parent-child relationship data structure
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP = convertToParentChildSSP(ssspTree);


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
    updateShortestPath(parentChildSSP, graphCSR, changedEdgesCSR, shortestDist);

    // Print the updated shortest path tree
    printShortestPathTree(parentChildSSP);

    return 0;
}
