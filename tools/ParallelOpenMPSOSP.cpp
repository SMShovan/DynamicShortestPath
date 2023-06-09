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

    // Read all lines first
    std::vector<std::string> lines;
    while (std::getline(inputFile, line)) {
        lines.push_back(line);
    }

    std::vector<std::vector<Edge>> csrMatrix(numRows);

    // Parallelize the processing of the lines
    #pragma omp parallel for
    for (int lineCount = 0; lineCount < lines.size(); ++lineCount) {
        std::istringstream lineStream(lines[lineCount]);
        int row, col;
        double value;
        if (!(lineStream >> row >> col >> value)) {
            #pragma omp critical
            {
                std::cerr << "Error: Invalid line format at line " << lineCount + 2 << " in the input file." << std::endl;
            }
            exit(1);
        }
        if (row < 1 || row > numRows || col < 1 || col > numCols) {
            #pragma omp critical
            {
                std::cerr << "Error: Invalid vertex indices at line " << lineCount + 2 << " in the input file." << std::endl;
            }
            exit(1);
        }
        #pragma omp critical
        {
            csrMatrix[row - 1].emplace_back(row - 1, col - 1, value);
        }
    }

    return csrMatrix;
}

#include <omp.h>
#include <vector>
#include <limits>

std::vector<std::vector<int>> dijkstra(const std::vector<std::vector<Edge>>& graphCSR, int sourceNode, std::vector<double>& shortestDist, std::vector<int>& parentList) {
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

        #pragma omp parallel for reduction(min:minDist)
        for (int j = 0; j < numNodes; ++j) {
            if (!visited[j] && shortestDist[j] < minDist) {
                #pragma omp critical
                {
                    if (shortestDist[j] < minDist) {
                        minDist = shortestDist[j];
                        minDistNode = j;
                    }
                }
            }
        }

        // Mark the minimum distance node as visited
        visited[minDistNode] = true;

        // Update the distances of the neighboring nodes
        for (const Edge& edge : graphCSR[minDistNode]) {
            int neighbor = edge.destination;
            double weight = edge.weight;
            if (!visited[neighbor] && shortestDist[minDistNode] + weight < shortestDist[neighbor]) {
                #pragma omp critical
                {
                    if (shortestDist[minDistNode] + weight < shortestDist[neighbor]) {
                        shortestDist[neighbor] = shortestDist[minDistNode] + weight;
                    }
                }
            }
        }
    }

    // Build the shortest path tree based on the shortest distances
    std::vector<std::vector<int>> ssspTree(numNodes);
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] != std::numeric_limits<double>::max()) {
            int parent = i + 1;
            for (const Edge& edge : graphCSR[i]) {
                int child = edge.destination + 1;
                if (shortestDist[child - 1] == shortestDist[i] + edge.weight) {
                    #pragma omp critical
                    {
                        ssspTree[parent - 1].push_back(child);
                        parentList[child] = parent;
                    }
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
    affectedNodes[node] = true;
    affectedDelNodes[node] = false;
    shortestDist[node] = std::numeric_limits<double>::infinity();

    for (int child : parentChildSSP[node].second) {
        markSubtreeAffected(parentChildSSP, shortestDist, affectedNodes, affectedDelNodes,  child);
    }
}



std::vector<std::pair<int, std::vector<int>>> convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree) {
    int n = ssspTree.size();
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP(n);

    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        parentChildSSP[i].first = i + 1;  // Set the node as the parent in the parent-child SSP structure
        parentChildSSP[i].second = ssspTree[i];  // Set the child nodes
    }

    // Optionally, For debug output serial, or remove it
    // std::cout << "Parent-Child SSP Tree:\n";
    // for (const auto& node : parentChildSSP) {
    //     std::cout << "Parent: " << node.first << ", Children: ";
    //     for (int child : node.second) {
    //         std::cout << child << " ";
    //     }
    //     std::cout << "\n";
    // }

    return parentChildSSP;
}




void updateShortestPath(std::vector<std::pair<int, std::vector<int>>>& ssspTree, std::vector<std::vector<Edge>>& graphCSR, const std::vector<std::vector<Edge>>& changedEdgesCSR, std::vector<double>& shortestDist, std::vector<int>& parentList) {
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
        std::cout << "Source: " << edge.source + 1 << " Destination: " << edge.destination + 1 << " Weight: " << edge.weight << "\n";
        
        int x,y; 
        // if(shortestDist[edge.source] > shortestDist[edge.destination]){
        //     x = edge.destination; 
        //     y = edge.source;
        // }  
        // else{
        //     x = edge.source; 
        //     y = edge.destination;
        // }
        x = edge.source;
        y = edge.destination;
        if (shortestDist[y ] > shortestDist[x] + edge.weight) {
            shortestDist[y] = shortestDist[x] + edge.weight;
            
            // Add the new edge to the graph
            graphCSR[x].push_back(Edge(x,y,edge.weight));

            int oldParent = parentList[y + 1];

            for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
            {
                
                if (ssspTree[oldParent - 1].second[j] == y + 1 )
                {
                    ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                    
                    break;
                }
            }
            parentList[y+1] = x + 1;
            ssspTree[x].second.push_back(y+1);
            affectedNodes[y] = true; 
            
            
        }
    }

    // Print the deleted edges and mark affected nodes
    std::cout << "Deleted Edges:\n";
    for (const Edge& edge : deletedEdges) {

        //Delete the edge from the graph 
        for (int i = 0; i < graphCSR[edge.source].size(); ++i) {
        if (graphCSR[edge.source][i].source == edge.source && graphCSR[edge.source][i].destination == edge.destination && graphCSR[edge.source][i].weight == edge.weight) {
            graphCSR[edge.source].erase(graphCSR[edge.source].begin() + i);
            break;
        }
    }

        std::cout << "Source: " << edge.source + 1 << " Destination: " << edge.destination + 1 << " Weight: " << edge.weight << "\n";
        affectedNodes[edge.destination] = true;
        affectedDelNodes[edge.destination] = true;
        shortestDist[edge.destination] = std::numeric_limits<double>::infinity();

        int oldParent = parentList[edge.destination + 1];

            for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
            {
                
                if (ssspTree[oldParent - 1].second[j] == edge.destination + 1 )
                {
                    ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                    
                    break;
                }
            }
            parentList[edge.destination + 1] =  -1;

        markSubtreeAffected(ssspTree, shortestDist, affectedNodes, affectedDelNodes, edge.destination);
    }

    // std::cout << "Affected Nodes:\n";
    // for (int i = 0; i < affectedNodes.size(); i++) {
    //     if (affectedNodes[i]) {
    //         std::cout << "Node " << i + 1 << "\n";
    //     }
    // }

    //printShortestPathTree(ssspTree);

    bool hasAffectedNodes = true;
    while (hasAffectedNodes) {
        hasAffectedNodes = false;

        // Check for affected nodes and update distances
        for (int v = 0; v < ssspTree.size(); ++v) {
            if (affectedNodes[v]) {
                
                affectedNodes[v] = false;  // Reset affected flag
                
                std::cout<< "Effected nodes source " << v+1;

                for (const Edge& edge : graphCSR[v] ) {
                    int n = edge.destination;
                    
                    std::cout<< " -> destination"<< n+1 <<"";

                    // for insertion
                    if (shortestDist[n] > shortestDist[v] + edge.weight) {
                        
                        shortestDist[n] = shortestDist[v] + edge.weight;
                        std::cout<< " :Distance Improved ";
                        
                        int oldParent = parentList[n + 1];
                        std::cout<< "old parent: " << oldParent;
                        for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
                        {
                            if (ssspTree[oldParent - 1].second[j] == n + 1 )
                            {
                                ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                                break;
                            }
                        }

                        parentList[n + 1] = v + 1; 
                        std::cout<< "New parent: " << v + 1 <<"\n";

                        ssspTree[v].second.push_back(n+1); 
                        affectedNodes[n] = true;  // Mark as affected
                        hasAffectedNodes = true;  // Set flag for next iteration

                    }
                    // for deletion
                    if ( shortestDist[v] == std::numeric_limits<double>::infinity() ){
                        shortestDist[n] = std::numeric_limits<double>::max();
                        std::cout<< " :Distance Infinity ";
                        
                        int oldParent = parentList[n + 1];
                        std::cout<< "old parent: " << oldParent;
                        for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
                        {
                            if (ssspTree[oldParent - 1].second[j] == n + 1 )
                            {
                                ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                                break;
                            }
                        }

                        parentList[n + 1] = -1; 
                        

                        
                        affectedNodes[n] = true;  // Mark as affected
                        hasAffectedNodes = true;
                    }
                }
            }
        }
    }
    int numNodes = ssspTree.size();
    std::vector<std::vector<int>> ssspTree2(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] != std::numeric_limits<double>::max()) {
            int parent = i + 1;
            for (const Edge& edge : graphCSR[i]) {
                int child = edge.destination + 1;
                if (shortestDist[child - 1] == shortestDist[i] + edge.weight) {
                    ssspTree2[parent - 1].push_back(child);
                }
            }
        }
    }
    // Print shortestDist
    std::cout << "\nShortest Distances:\n";
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] == std::numeric_limits<double>::max()) {
            std::cout << "Node " << i + 1 << ": Infinity\n";
        } else {
            std::cout << "Node " << i + 1 << ": " << shortestDist[i] << "\n";
        }
    }

    // Print ssspTree
    std::cout << "Correct Shortest Path Tree:\n";
    for (int i = 0; i < numNodes; ++i) {
        std::cout << "Node " << i + 1 << ": ";
        for (int child : ssspTree2[i]) {
            std::cout << child << " ";
        }
        std::cout << "\n";
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
    std::vector<double> shortestDist(graphCSR.size(), std::numeric_limits<double>::max());
    std::vector<int> parent(graphCSR.size() + 1, -1);
    std::vector<std::vector<int>> ssspTree = dijkstra(graphCSR, sourceNode, shortestDist, parent);

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
    updateShortestPath(parentChildSSP, graphCSR, changedEdgesCSR, shortestDist, parent);

    // Print the updated shortest path tree
    printShortestPathTree(parentChildSSP);

    return 0;
}
