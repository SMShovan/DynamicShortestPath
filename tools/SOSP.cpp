#include <iostream>
#include <fstream>
#include <vector>
#include "csr_conversion.h"
#include "changedEdgesToCSR.h"
#include "dijkstra.h"

// Function to update the shortest path based on the changed edges
void updateShortestPath(std::vector<double>& shortestPath, const std::vector<std::vector<Edge>>& changedEdgesCSR) {
    std::vector<bool> affectedVertices(shortestPath.size(), false);

    // Step 1: Identify affected vertices
    for (const std::vector<Edge>& edges : changedEdgesCSR) {
        for (const Edge& edge : edges) {
            int source = edge.source;
            int destination = edge.destination;

            affectedVertices[source] = true;
            affectedVertices[destination] = true;
        }
    }

    // Step 2: Update CSR matrix and shortest path tree
    for (const std::vector<Edge>& edges : changedEdgesCSR) {
        for (const Edge& edge : edges) {
            int source = edge.source;
            int destination = edge.destination;
            double weight = edge.weight;

            // Check if both vertices are affected
            if (affectedVertices[source] && affectedVertices[destination]) {
                // Update CSR matrix
                shortestPath[destination] = weight;
            }
            // Check if only the source vertex is affected
            else if (affectedVertices[source]) {
                // Update CSR matrix
                for (Edge& e : graphCSR[source]) {
                    if (e.destination == destination) {
                        e.weight = weight;
                        break;
                    }
                }

                // Update shortest path tree (add or remove edge)
                if (weight > 0) {
                    // Add edge to the shortest path tree if it doesn't exist
                    bool edgeExists = false;
                    for (const Edge& e : graphCSR[source]) {
                        if (e.destination == destination) {
                            edgeExists = true;
                            break;
                        }
                    }
                    if (!edgeExists) {
                        graphCSR[source].push_back(Edge(destination, weight));
                    }
                } else {
                    // Remove edge from the shortest path tree if it exists
                    for (auto it = graphCSR[source].begin(); it != graphCSR[source].end(); ++it) {
                        if (it->destination == destination) {
                            graphCSR[source].erase(it);
                            break;
                        }
                    }
                }
            }
            // Check if only the destination vertex is affected
            else if (affectedVertices[destination]) {
                // Update CSR matrix
                for (Edge& e : graphCSR[destination]) {
                    if (e.destination == source) {
                        e.weight = weight;
                        break;
                    }
                }

                // Update shortest path tree (add or remove edge)
                if (weight > 0) {
                    // Add edge to the shortest path tree if it doesn't exist
                    bool edgeExists = false;
                    for (const Edge& e : graphCSR[destination]) {
                        if (e.destination == source) {
                            edgeExists = true;
                            break;
                        }
                    }
                    if (!edgeExists) {
                        graphCSR[destination].push_back(Edge(source, weight));
                    }
                } else {
                    // Remove edge from the shortest path tree if it exists
                    for (auto it = graphCSR[destination].begin(); it != graphCSR[destination].end(); ++it) {
                        if (it->destination == source) {
                            graphCSR[destination].erase(it);
                            break;
                        }
                    }
                }
            }
            // Neither vertex is affected, no changes required
        }
    }
}


int main(int argc, char** argv) {
    std::string graphFile;
    std::string changedEdgesFile;
    int sourceNode;

    // Parse command-line options
    // ...

    // Read the input graph in Matrix Market format
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
