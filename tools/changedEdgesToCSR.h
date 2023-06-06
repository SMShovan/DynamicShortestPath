#ifndef CHANGED_EDGES_TO_CSR_H
#define CHANGED_EDGES_TO_CSR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// Struct to represent a modified edge
struct ModifiedEdge {
    int source;
    int destination;
    double weight;
    int changeType; // +1 for insertion, -1 for deletion
};

// Function to convert modified edges to CSR format
bool convertToCSR(const std::string& inputFile, std::vector<int>& rowPointers, std::vector<int>& columnIndices, std::vector<double>& values) {
    std::ifstream file(inputFile);
    std::string line;

    if (!file.is_open()) {
        std::cout << "Unable to open input file: " << inputFile << std::endl;
        return false;
    }

    // Skip comment lines
    while (std::getline(file, line) && line[0] == '%') {}

    // Read matrix dimensions and number of modified edges
    int numRows, numCols, numModifiedEdges;
    std::istringstream dimensions(line);
    dimensions >> numRows >> numCols >> numModifiedEdges;

    rowPointers.resize(numRows + 1, 0);
    columnIndices.reserve(numModifiedEdges);
    values.reserve(numModifiedEdges);

    int numEntries = 0;

    // Read modified edges and construct CSR matrix
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        ModifiedEdge edge;
        iss >> edge.source >> edge.destination >> edge.weight >> edge.changeType;

        rowPointers[edge.source]++;
        columnIndices.push_back(edge.destination);
        values.push_back(edge.weight);
        numEntries++;
    }

    // Update row pointers to cumulative sum
    for (int i = 1; i <= numRows; i++) {
        rowPointers[i] += rowPointers[i - 1];
    }

    file.close();

    return true;
}

#endif  // CHANGED_EDGES_TO_CSR_H

