#ifndef CHANGED_EDGES_TO_CSR_H
#define CHANGED_EDGES_TO_CSR_H

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

std::vector<std::vector<Edge>> convertToCSR(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open input file: " << filename << std::endl;
        return std::vector<std::vector<Edge>>();
    }

    // Read the Matrix Market file and convert it to CSR format
    std::string line;
    int numRows, numCols, numNonZero;
    std::getline(inputFile, line);
    std::istringstream headerStream(line);
    headerStream >> numRows >> numCols >> numNonZero;

    std::vector<std::vector<Edge>> csrMatrix(numRows);

    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        int row, col, change;
        double value;
        lineStream >> row >> col >> value >> change;
        csrMatrix[row - 1].emplace_back(row - 1, col - 1, value);
    }

    inputFile.close();

    return csrMatrix;
}

#endif  // CHANGED_EDGES_TO_CSR_H
