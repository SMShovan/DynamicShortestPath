#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>
#include <chrono>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cmath>
#include <omp.h>
int NUM_THREADS = 4;

using namespace std;

template<typename T>
void removeElementsFromQueue(std::queue<T>& q, const T& value) {
    std::queue<T> tempQueue;

    while (!q.empty()) {
        T frontElement = q.front();
        q.pop();

        // Skip the elements with the specified value
        if (frontElement != value) {
            tempQueue.push(frontElement);
        }
    }

    // Swap the contents back to the original queue
    std::swap(q, tempQueue);
}

struct Edge {
    int source;
    int destination;
    double weight;

    Edge(int src, int dest, double w) : source(src), destination(dest), weight(w) {}

    // operator < needed to keep edges in the set
    // bool operator<(const Edge& other) const {
    //     return this->weight < other.weight;
    // }
};

double find75Percentile(std::vector<int> &data) {
    if (data.empty()) {
        throw std::invalid_argument("Cannot find percentile of an empty vector.");
    }

    // Sort the vector
    std::sort(data.begin(), data.end());

    // Calculate index (75th percentile)
    double idx = 0.75 * (data.size() - 1);

    // Calculate the value at the 75th percentile position
    // Note: We need to account for fractional indices by doing linear interpolation.
    if (std::floor(idx) == idx) { // If the index is an integer
        return data[static_cast<int>(idx)];
    } else { // If the index is a fraction
        double frac = idx - std::floor(idx); // Fractional part of the index
        return (1.0 - frac) * data[static_cast<int>(idx)] + frac * data[static_cast<int>(idx) + 1];
    }
}

std::vector<std::vector<Edge>> predecessor;
std::vector<int> predecessorCount;
bool operator==(const Edge& lhs, const Edge& rhs) {
    return lhs.source == rhs.source && lhs.destination == rhs.destination && lhs.weight == rhs.weight;
}

std::vector<std::vector<Edge>> convertToCSR(std::ifstream& inputFile, bool isGraph) {
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
        if(isGraph)
        {
            predecessor.resize(numCols);
            predecessorCount.resize(numCols, 0);
            predecessor[col - 1].emplace_back(row, col, value);
            predecessorCount[col - 1]++;
        }
            
    }

    // std::cout << "predecessorCount: ";
    // for (int val : predecessorCount) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    // if (isGraph)
    // {
    //     std::cout << "In-Degree Matrix:" << std::endl;
    //     for (int col = 0; col < predecessor.size(); ++col) {
    //         for (const auto& edge : predecessor[col]) {
    //             std::cout << "(" << edge.source << ", " << edge.destination << ", " << edge.weight << ") ";
    //         }
    //         std::cout << std::endl;
    //     }  
    // }    
        



    return csrMatrix;
}

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
        double minDist = std::numeric_limits<double>::infinity();
        for (int j = 0; j < numNodes; ++j) {
            if (!visited[j] && shortestDist[j] < minDist) {
                minDist = shortestDist[j];
                minDistNode = j;
            }
        }

        // Mark the minimum distance node as visited
        if (minDistNode == -1) {
            break;
        }
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
    std::vector<bool> cycleCheck (numNodes, false);
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] != std::numeric_limits<double>::infinity()) {
            int parent = i + 1;
            for (const Edge& edge : graphCSR[i]) {
                int child = edge.destination + 1;
                if (shortestDist[child - 1] == shortestDist[i] + edge.weight && !cycleCheck[child - 1]) {
                    ssspTree[parent - 1].push_back(child);
                    cycleCheck[child - 1 ] = true;
                    parentList[child] = parent;
                }
            }
        }
    }

    // // Print shortestDist
    // std::cout << "Parent list:\n";
    // for (int i = 1; i <= numNodes; ++i) {
    //     std::cout<< "Child: "<< i << " Parent: "<< parentList[i];
    // }

    // // Print shortestDist
    // std::cout << "Shortest Distances:\n";
    // for (int i = 0; i < numNodes; ++i) {
    //     if (shortestDist[i] == std::numeric_limits<double>::infinity()) {
    //         std::cout << "Node " << i + 1 << ": Infinity\n";
    //     } else {
    //         std::cout << "Node " << i + 1 << ": " << shortestDist[i] << "\n";
    //     }
    // }

    // // Print ssspTree
    // std::cout << "Shortest Path Tree:\n";
    // for (int i = 0; i < numNodes; ++i) {
    //     std::cout << "Node " << i + 1 << ": ";
    //     for (int child : ssspTree[i]) {
    //         std::cout << child << " ";
    //     }
    //     std::cout << "\n";
    // }

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


void markSubtreeAffected(std::vector<std::pair<int, std::vector<int>>>& parentChildSSP, std::vector<double>& shortestDist, std::vector<bool>& affectedNodes, std::queue<int>& affectedNodesQueue, std::vector<bool>& affectedDelNodes, int node) {
    affectedNodes[node] = true;
    affectedNodesQueue.push(node);
    affectedDelNodes[node] = true;
    shortestDist[node] = std::numeric_limits<double>::infinity();

    #pragma omp parallel for
    for (size_t i = 0; i < parentChildSSP[node].second.size(); ++i) {
        int child = parentChildSSP[node].second[i];
        if (!affectedDelNodes[child]) {
            markSubtreeAffected(parentChildSSP, shortestDist, affectedNodes, affectedNodesQueue, affectedDelNodes, child);
        }
    }
}

std::vector<std::pair<int, std::vector<int>>> convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree) {
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP(ssspTree.size());

    for (int i = 0; i < ssspTree.size(); ++i) {
        parentChildSSP[i].first = i + 1;  // Set the node as the parent in the parent-child SSP structure
        parentChildSSP[i].second = ssspTree[i];  // Set the child nodes
    }

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
    std::queue<int> affectedNodesQueue;
    int threshold;

    std::unordered_map<int, std::vector<Edge>> insertedEdgeMap;
    std::unordered_map<int, std::vector<Edge>> deletedEdgeMap;

    // Identify inserted and deleted edges
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (size_t i = 0; i < changedEdgesCSR.size(); ++i) {
        std::vector<Edge> privateInsertedEdges;
        std::vector<Edge> privateDeletedEdges;


        const std::vector<Edge>& row = changedEdgesCSR[i];
        for (size_t j = 0; j < row.size(); ++j) {
            const Edge& edge = row[j];
            if (edge.weight > 0) {
                privateInsertedEdges.push_back(edge);
            } else if (edge.weight < 0) {
                privateDeletedEdges.push_back(edge);
            }
        }

        #pragma omp critical
        {
            // Merge privateInsertedEdges into insertedEdges
            insertedEdges.insert(insertedEdges.end(), privateInsertedEdges.begin(), privateInsertedEdges.end());

            // Merge privateDeletedEdges into deletedEdges
            deletedEdges.insert(deletedEdges.end(), privateDeletedEdges.begin(), privateDeletedEdges.end());
        }
    }

        

    // for (const auto& pair : insertedEdgeMap) {
    //     std::cout << "Edges ending at node " << pair.first << ":\n";
    //     for (const auto& edge : pair.second) {
    //         std::cout << "Start: " << edge.source << ", End: " << edge.destination << ", Weight: " << edge.weight << "\n";
    //     }
    // }

    

    // Print the inserted edges and identify affected nodes
    //std::cout << "Inserted Edges:\n";
    
    // #pragma omp parallel for
    // for (int i = 0; i < insertedEdgeMap.size(); i++) 
    //     for (const Edge& edge : insertedEdgeMap[i]) {
    //         //std::cout << "Source: " << edge.source + 1 << " Destination: " << edge.destination + 1 << " Weight: " << edge.weight << "\n";
            
    //         int x,y; 
    //         // if(shortestDist[edge.source] > shortestDist[edge.destination]){
    //         //     x = edge.destination; 
    //         //     y = edge.source;
    //         // }  
    //         // else{
    //         //     x = edge.source; 
    //         //     y = edge.destination;
    //         // }
    //         x = edge.source;
    //         y = edge.destination;
    //         if (shortestDist[y] > shortestDist[x] + edge.weight) {
    //             shortestDist[y] = shortestDist[x] + edge.weight;
                
    //             // Add the new edge to the graph
    //             graphCSR[x].push_back(Edge(x,y,edge.weight));

    //             int oldParent = parentList[y + 1];

    //             for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
    //             {
                    
    //                 if (ssspTree[oldParent - 1].second[j] == y + 1 )
    //                 {
    //                     ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                        
    //                     break;
    //                 }
    //             }
    //             parentList[y+1] = x + 1;
    //             ssspTree[x].second.push_back(y+1);
    //             affectedNodes[y] = true; 
    //             affectedNodesQueue.push(y);
                
                
    //         }
    //     }

    #pragma omp parallel for 
    for (int i = 0; i < insertedEdgeMap.size(); i++) {
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < insertedEdgeMap[i].size(); j++) {
            const Edge& edge = insertedEdgeMap[i][j];

            int x, y;
            x = edge.source;
            y = edge.destination;

            int newShortestDist = shortestDist[x] + edge.weight;

            if (newShortestDist < shortestDist[y]) {
                // Using atomic operation to update shortestDist[y] without a critical section
                #pragma omp atomic write
                shortestDist[y] = newShortestDist;

                // Add the new edge to the graph (no critical section needed)
                #pragma omp critical
                {
                graphCSR[x].push_back(Edge(x, y, edge.weight));
                }

                int oldParent = parentList[y + 1];

                // Task-level parallelism for erasing ssspTree elements
                #pragma omp task
                {
                    auto it = std::find(ssspTree[oldParent - 1].second.begin(), ssspTree[oldParent - 1].second.end(), y + 1);
                    if (it != ssspTree[oldParent - 1].second.end()) {
                        // Using atomic operation to erase ssspTree element without a critical section
                        #pragma omp atomic write
                        *it = -1; // Replace the value with a special marker value (-1) for later removal
                    }
                }

                parentList[y + 1] = x + 1;

                // Task-level parallelism for adding ssspTree elements
                #pragma omp task
                {
                    #pragma omp critical
                    ssspTree[x].second.push_back(y + 1);
                }

                // Task-level parallelism for updating affectedNodes and affectedNodesQueue
                #pragma omp task
                {
                    #pragma omp critical
                    {
                        affectedNodes[y] = true;
                        affectedNodesQueue.push(y);
                    }
                }
            }
        }
    }
    #pragma omp taskwait


    

    // Print the deleted edges and mark affected nodes
    //std::cout << "Deleted Edges:\n";

    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < deletedEdgeMap.size(); i++)
        for (const Edge& edge : deletedEdges) {

            //Delete the edge from the graph 
            for (int i = 0; i < graphCSR[edge.source].size(); ++i) {
            if (graphCSR[edge.source][i].source == edge.source && graphCSR[edge.source][i].destination == edge.destination) {
                graphCSR[edge.source].erase(graphCSR[edge.source].begin() + i);
                break;
            }
            }

            // After deletion from graph, check if the deleted edges belong to ssspTree. 

            bool inTree = false;

        // Delete the element from the predecessor as well

            for (auto& preEdge : predecessor[edge.destination]) 
            if (edge.source + 1 == preEdge.source) {
                
                // Find the iterator pointing to the element
                auto it = std::find(predecessor[edge.destination].begin(), predecessor[edge.destination].end(), preEdge);
                
                for (int i = 0; i < ssspTree[edge.source].second.size() ; i++)
                {
                    if (ssspTree[edge.source].second[i] == edge.destination + 1)
                    //std::cout<<"Found in tree"<<std::endl;
                    inTree = true;

                }

                if (it != predecessor[edge.destination].end()) {

                    predecessor[edge.destination].erase(it);
                    predecessorCount[edge.destination]--;
                    
                }
            }

        threshold = find75Percentile(predecessorCount);
        

        if (!inTree)
            continue;

        

        

        markSubtreeAffected(ssspTree, shortestDist, affectedNodes, affectedNodesQueue, affectedDelNodes, edge.destination);

        

        //std::cout << "Source: " << edge.source + 1 << " Destination: " << edge.destination + 1 << " Weight: " << edge.weight << "\n";
        affectedNodes[edge.destination] = true;
        affectedNodesQueue.push(edge.destination);
        affectedDelNodes[edge.destination] = true;
        shortestDist[edge.destination] = std::numeric_limits<double>::infinity();

        // find new parent or put the shortest distance to infinity
        int n = edge.destination; 
        //std::cout<<"\nPredecessor of "<< n + 1 <<" ";
        int newDistance = std::numeric_limits<double>::infinity();
        int newParentIndex = -1; 


        

        for (int i = 0; i < predecessor[n].size(); i++)
        {
            if (shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight < newDistance)
            {
                newDistance = shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight;

                newParentIndex = predecessor[n][i].source - 1;
            } 
        }
        int oldParent = parentList[edge.destination + 1];
        

        if (newParentIndex == -1)
        {
            parentList[n + 1] = -1; 
            shortestDist[n] = std::numeric_limits<double>::infinity();

        }   
        else 
        {              
            ssspTree[newParentIndex].first = newParentIndex + 1; 
            shortestDist[n] = newDistance; 
            //std::cout<< "New parent: "<< newParentIndex + 1<< " child "<< n + 1<< "\n"; 
            ssspTree[newParentIndex].second.push_back(n+1);      
        }

    

        //std::cout << "Old parent" <<  oldParent << "\n";

            for (int j = 0; j < ssspTree[oldParent - 1].second.size() && oldParent != -1 ; j++ )
            {
                
                if (ssspTree[oldParent - 1].second[j] == edge.destination + 1 )
                {
                    ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                    break;
                }
            }

        
        parentList[edge.destination + 1] = newParentIndex + 1; 

        
            

    }

     

    // std::cout << "Affected Nodes:\n";
    // for (int i = 0; i < affectedNodes.size(); i++) {
    //     if (affectedNodes[i]) {
    //         std::cout << "Node " << i + 1 << "\n";
    //     }
    // }

    //printShortestPathTree(ssspTree);
    
#pragma omp parallel num_threads(NUM_THREADS)
{

    bool hasAffectedNodes = true;
    while (hasAffectedNodes) {


        hasAffectedNodes = false;

        // Check for affected nodes and update distances
        while(!affectedNodesQueue.empty()){
        //for (int v = 0; v < ssspTree.size(); ++v) {
            int v = affectedNodesQueue.front();
            affectedNodesQueue.pop();
            if (affectedNodes[v]) {

                affectedNodes[v] = false;  // Reset affected flag
                

                
                //std::cout<< "Effected nodes source " << v+1;
                #pragma omp parallel for if(predecessorCount[v] <= threshold)
                for (size_t i = 0; i < graphCSR[v].size(); ++i) {
                    const Edge& edge = graphCSR[v][i];
                    int n = edge.destination;

                    // for insertion
                    if (shortestDist[n] > shortestDist[v] + edge.weight && shortestDist[v] != std::numeric_limits<double>::infinity()) {
                        shortestDist[n] = shortestDist[v] + edge.weight;

                        int oldParent = parentList[n + 1];
                        #pragma omp parallel for if(predecessorCount[v] > threshold)
                        for (size_t j = 0; j < ssspTree[oldParent - 1].second.size(); j++) {
                            if (ssspTree[oldParent - 1].second[j] == n + 1) {
                                ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                            }
                        }

                        parentList[n + 1] = v + 1;
                        ssspTree[v].second.push_back(n + 1);
                        affectedNodes[n] = true;  // Mark as affected
                        affectedNodesQueue.push(n);
                        hasAffectedNodes = true;  // Set flag for next iteration
                    }

                    // for deletion
                    if (shortestDist[v] == std::numeric_limits<double>::infinity()) {
                        int newDistance = std::numeric_limits<double>::infinity();
                        int newParentIndex = -1;
                        #pragma omp parallel for if(predecessorCount[v] > threshold)
                        for (size_t i = 0; i < predecessor[n].size(); i++) {
                            if (shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight < newDistance) {
                                newDistance = shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight;
                                newParentIndex = predecessor[n][i].source - 1;
                            }
                        }
                        if (n + 1 == 1)
                            continue;
                        int oldParent = parentList[n + 1];

                        if (newParentIndex == -1) {
                            parentList[n + 1] = -1;
                            shortestDist[n] = std::numeric_limits<double>::infinity();
                        } else {
                            ssspTree[newParentIndex].first = newParentIndex + 1;
                            shortestDist[n] = newDistance;
                            ssspTree[newParentIndex].second.push_back(n + 1);
                        }

                        #pragma omp parallel for if(predecessorCount[v] > threshold)
                        for (size_t j = 0; j < ssspTree[oldParent - 1].second.size(); j++) {
                            if (ssspTree[oldParent - 1].second[j] == n + 1) {
                                ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                            }
                        }
                        parentList[n + 1] = newParentIndex + 1;
                        affectedNodes[n] = true;  // Mark as affected
                        affectedNodesQueue.push(n);
                        hasAffectedNodes = true;
                    }
                }
            }
        }
    }
}

    


    

    // int numNodes = ssspTree.size();
    // std::vector<std::vector<int>> ssspTree2(numNodes);
    // std::vector<bool> cycleCheck(numNodes, false);
    // for (int i = 0; i < numNodes; ++i) {
    //     if (shortestDist[i] != std::numeric_limits<double>::infinity()) {
    //         int parent = i + 1;
    //         for (const Edge& edge : graphCSR[i]) {
    //             int child = edge.destination + 1;
    //             if (shortestDist[child - 1] == shortestDist[i] + edge.weight && !cycleCheck[child - 1]) {
    //                 ssspTree2[parent - 1].push_back(child);
    //                 cycleCheck[child - 1] = true; 
    //             }
    //         }
    //     }
    // }
    // // Print shortestDist
    // std::cout << "\nShortest Distances:\n";
    // for (int i = 0; i < numNodes; ++i) {
    //     if (shortestDist[i] == std::numeric_limits<double>::infinity()) {
    //         std::cout << "Node " << i + 1 << ": Infinity\n";
    //     } else {
    //         std::cout << "Node " << i + 1 << ": " << shortestDist[i] << "\n";
    //     }
    // }

    // //Print ssspTree
    // std::cout << "Correct Shortest Path Tree:\n";
    // for (int i = 0; i < numNodes; ++i) {
    //     std::cout << "Node " << i + 1 << ": ";
    //     for (int child : ssspTree2[i]) {
    //         std::cout << child << " ";
    //     }
    //     std::cout << "\n";
    // }

}

std::vector<std::vector<std::pair<int, int>>>  combineGraphs(const std::vector<std::pair<int, std::vector<int>>>& parentChildSSP1, const std::vector<std::pair<int, std::vector<int>>>& parentChildSSP2) {

     int maxNodeIndex = 0;

    // Find the maximum node index in both input graphs
    for (const auto& node : parentChildSSP1) {
        maxNodeIndex = std::max(maxNodeIndex, node.first);
        for (int child : node.second) {
            maxNodeIndex = std::max(maxNodeIndex, child);
        }
    }

    

    for (const auto& node : parentChildSSP2) {
        maxNodeIndex = std::max(maxNodeIndex, node.first);
        for (int child : node.second) {
            maxNodeIndex = std::max(maxNodeIndex, child);
        }
    }

    

    // Initialize the combinedGraph vector with the appropriate size
    std::vector<std::vector<std::pair<int, int>>> combinedGraph(maxNodeIndex);


    // Iterate over nodes in parentChildSSP1
    for (const auto& node : parentChildSSP1) {
        int parent1 = node.first;
        const std::vector<int>& children1 = node.second;

        // Add children of parent1 to combinedGraph with edge weight 1
        for (int child1 : children1) {
            combinedGraph[parent1 - 1].push_back({child1,2});
        }
    }

    
    
    // Iterate over nodes in parentChildSSP2
    for (const auto& node : parentChildSSP2) {
        int parent2 = node.first;
        const std::vector<int>& children2 = node.second;

        for (int child2 : children2) {
            // Check if the child node already exists in combinedGraph
            bool foundChild = false;
            for (auto& childPair : combinedGraph[parent2 - 1]) {
                if (childPair.first == child2) {
                    foundChild = true;
                    childPair.second = 1;
                    break;
                }
            }

            // Add the child node with the appropriate edge weight
            if (!foundChild) {
                combinedGraph[parent2 - 1].push_back({child2, 2});
            }
        }
    }
    
    // Print the combined graph
    // std::cout << "\nCombined Graph:\n";
    // for (int i = 0; i < combinedGraph.size(); i++) {
    //     std::cout << "Node " << i + 1 << ": ";
    //     for ( auto& childPair : combinedGraph[i]) {
    //         int child = childPair.first;
    //         std::cout << "(" << child << ", " << childPair.second << ") ";
    //     }
    //     std::cout << "\n";
    // }

    return combinedGraph;
}

std::vector<std::pair<int, int>> dijkstra2(const std::vector<std::vector<std::pair<int, int>>>& graph, int source) {
    int numNodes = graph.size();
    std::vector<int> distance(numNodes, std::numeric_limits<int>::max());
    std::vector<int> parent(numNodes, -1);
    std::vector<bool> visited(numNodes, false);

    distance[source - 1] = 0;

    // Use a priority queue to get the node with the smallest distance
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
    pq.push({0, source - 1});

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (visited[u]) {
            continue;
        }

        visited[u] = true;

        for (const auto& edge : graph[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (!visited[v - 1] && distance[u] + weight < distance[v - 1]) {
                distance[v - 1] = distance[u] + weight;
                parent[v - 1] = u;
                pq.push({distance[v - 1], v - 1});
            }
        }
    }

    std::vector<std::pair<int, int>> shortestPathTree;
    for (int i = 0; i < numNodes; ++i) {
        if (i != source - 1) {
            shortestPathTree.emplace_back(i + 1, parent[i] + 1);
        }
    }

    return shortestPathTree;
}

std::vector<std::pair<int, int>> bellmanFord(const std::vector<std::vector<std::pair<int, int>>>& combinedGraph, int source = 1) {
    int numNodes = combinedGraph.size();
    std::vector<int> distance(numNodes, std::numeric_limits<int>::max());
    std::vector<int> parent(numNodes, -1);

    distance[source - 1] = 0;

    

    // Relax edges repeatedly
    for (int i = 1; i < numNodes; ++i) {
        for (int u = 0; u < numNodes; ++u) {
            
            for (const auto& edge : combinedGraph[u]) {
                std::cout<<"Good till here " <<std::endl;
                
                int v = edge.first;
                int weight = edge.second;
                    
                if (distance[u] != std::numeric_limits<int>::max() && distance[u] + weight < distance[v - 1]) {
                    distance[v - 1] = distance[u] + weight;
                    parent[v - 1] = u;
                }
            }
        }
    }
    

    // //Check for negative cycles
    // for (int u = 0; u < numNodes; ++u) {
    //     for (const auto& edge : combinedGraph[u]) {
    //         int v = edge.first;
    //         int weight = edge.second;
    //         if (distance[u] != std::numeric_limits<int>::max() && distance[u] + weight < distance[v - 1]) {
    //             std::cout << "Graph contains a negative cycle!" << std::endl;
                
    //         }
    //     }
    //}

    std::vector<std::pair<int, int>> shortestPathTree;
    for (int i = 0; i < numNodes; ++i) {
        if (i + 1 != source) { // Exclude the source node
            shortestPathTree.emplace_back(i + 1, parent[i] + 1);
        }
    }
    return shortestPathTree;
}

std::vector<std::pair<int, std::vector<int>>> convertParentCombinedToParentChildSSP(const std::vector<int>& parentCombined) {
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP(parentCombined.size());

    for (int i = 0; i < parentCombined.size(); ++i) {
        parentChildSSP[i].first = i + 1;  // Set the node as the parent in the parent-child SSP structure
        int parent = parentCombined[i];
        if (parent != -1) {
            parentChildSSP[parent].second.push_back(i + 1);  // Set the child nodes
        }
    }

    return parentChildSSP;
}


int main(int argc, char** argv) {
    //std::cout<< "Parallel code started"<< std::endl;
    // Check the command-line arguments
    if (argc < 10) {
        std::cerr << "Usage: ./program -g <graph_file> -c <changed_edges_file> -s <source_node>\n";
        return 1;
    }

    std::string graphFile1;
    std::string changedEdgesFile1;
    std::string graphFile2;
    std::string changedEdgesFile2;
    int sourceNode;
    

    // Parse the command-line arguments
    for (int i = 1; i < argc; i += 2) {
        std::string option(argv[i]);
        std::string argument(argv[i + 1]);

        if (option == "-g1") {
            graphFile1 = argument;
        } else if (option == "-c1") {
            changedEdgesFile1 = argument;
        } else if (option == "-s") {
            sourceNode = std::stoi(argument);
        }else if (option == "-t") {
            NUM_THREADS = std::stoi(argument);
        }
        else if (option == "-g2") {
            graphFile2 = argument;
        }
        else if (option == "-c2") {
            changedEdgesFile2 = argument;
        } else {
            std::cerr << "Invalid option: " << option << "\n";
            return 1;
        }
    }

    std::ifstream graphInputFile1(graphFile1);
    if (!graphInputFile1.is_open()) {
        std::cerr << "Unable to open graph file 1: " << graphFile1 << std::endl;
        return 1;
    }

    std::ifstream graphInputFile2(graphFile2);
    if (!graphInputFile2.is_open()) {
        std::cerr << "Unable to open graph file 2: " << graphFile2 << std::endl;
        return 1;
    }

    // Convert the input graphs to CSR format
    std::cout<<"Graph is converting to CSR"<<std::endl;
    std::vector<std::vector<Edge>> graphCSR1 = convertToCSR(graphInputFile1, true);
    std::vector<std::vector<Edge>> graphCSR2 = convertToCSR(graphInputFile2, true);
    graphInputFile1.close();
    graphInputFile2.close();

    // Find the initial shortest path trees using Dijkstra's algorithm
    std::cout<<"Finding the shortest distance of graph 1"<<std::endl;
    std::vector<double> shortestDist1(graphCSR1.size(), std::numeric_limits<double>::infinity());
    std::vector<int> parent1(graphCSR1.size() + 1, -1);
    std::cout<<"Dijkstra started for of graph 1"<<std::endl;
    std::vector<std::vector<int>> ssspTree1 = dijkstra(graphCSR1, sourceNode, shortestDist1, parent1);
    std::cout<<"Finding the shortest distance of graph 2"<<std::endl;
    std::vector<double> shortestDist2(graphCSR2.size(), std::numeric_limits<double>::infinity());
    std::vector<int> parent2(graphCSR2.size() + 1, -1);
    std::cout<<"Dijkstra started for of graph 1"<<std::endl;
    std::vector<std::vector<int>> ssspTree2 = dijkstra(graphCSR2, sourceNode, shortestDist2, parent2);

    // Convert the ssspTrees to parent-child relationship data structure
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP1 = convertToParentChildSSP(ssspTree1);
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP2 = convertToParentChildSSP(ssspTree2);

    std::ifstream changedEdgesInputFile1(changedEdgesFile1);
    if (!changedEdgesInputFile1.is_open()) {
        std::cerr << "Unable to open changed edges file 1: " << changedEdgesFile1 << std::endl;
        return 1;
    }

    std::ifstream changedEdgesInputFile2(changedEdgesFile2);
    if (!changedEdgesInputFile2.is_open()) {
        std::cerr << "Unable to open changed edges file 2: " << changedEdgesFile2 << std::endl;
        return 1;
    }

    // Convert the changed edges to CSR format
    std::vector<std::vector<Edge>> changedEdgesCSR1 = convertToCSR(changedEdgesInputFile1, false);
    std::vector<std::vector<Edge>> changedEdgesCSR2 = convertToCSR(changedEdgesInputFile2, false);
    changedEdgesInputFile1.close();
    changedEdgesInputFile2.close();

    // Update the shortest path trees based on the changed edges
    std::cout<<"Updating the shortest path for of graph 1"<<std::endl;
    updateShortestPath(parentChildSSP1, graphCSR1, changedEdgesCSR1, shortestDist1, parent1);
    std::cout<<"Updating the shortest path for of graph 2"<<std::endl;
    updateShortestPath(parentChildSSP2, graphCSR2, changedEdgesCSR2, shortestDist2, parent2);

    // Print the updated shortest path trees
    //std::cout << "Shortest Path Tree 1:\n";
    printShortestPathTree(parentChildSSP1);

    //std::cout << "\nShortest Path Tree 2:\n";
    printShortestPathTree(parentChildSSP2);

    std::cout<<"Combining the graphs"<<std::endl;
    std::vector<std::vector<std::pair<int, int>>>  combinedGraph = combineGraphs(parentChildSSP1, parentChildSSP2);

    std::cout<<"Dijkstra is running"<<std::endl;

    std::vector<std::pair<int, int>> shortestPathTree = dijkstra2(combinedGraph, 1);

    //Print the shortest path tree
    for (const auto& node : shortestPathTree) {
        std::cout << "Node " << node.first << ": Parent = " << node.second << "\n";
    }



    std::cout<<"DONE"<<std::endl;

    return 0;
}
