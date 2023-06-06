#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <getopt.h>

#define INF std::numeric_limits<double>::infinity()

// Edge struct to represent an edge in the graph
struct Edge {
    int destination;
    double weight;
};

// Function to find the shortest path using Dijkstra's algorithm
std::vector<double> dijkstra(const std::vector<std::vector<Edge>>& graph, int source) {
    int numNodes = graph.size();
    std::vector<double> distances(numNodes, INF);
    distances[source] = 0.0;

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    pq.push(std::make_pair(0.0, source));

    while (!pq.empty()) {
        int currentNode = pq.top().second;
        double currentDistance = pq.top().first;
        pq.pop();

        if (currentDistance > distances[currentNode]) {
            continue;
        }

        for (const Edge& edge : graph[currentNode]) {
            int nextNode = edge.destination;
            double weight = edge.weight;

            if (currentDistance + weight < distances[nextNode]) {
                distances[nextNode] = currentDistance + weight;
                pq.push(std::make_pair(distances[nextNode], nextNode));
            }
        }
    }

    return distances;
}

// Function to print the CSR matrix
void printCSR(const std::vector<std::vector<Edge>>& csrMatrix) {
    std::cout << "CSR Matrix:" << std::endl;
    for (int i = 0; i < csrMatrix.size(); i++) {
        for (const Edge& edge : csrMatrix[i]) {
            std::cout << i + 1 << " " << edge.destination + 1 << " " << edge.weight << std::endl;
        }
    }
}

#endif  // DIJKSTRA_H

