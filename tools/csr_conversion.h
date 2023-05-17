#ifndef CSR_CONVERSION_H
#define CSR_CONVERSION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

struct CSRMatrix {
    std::vector<int> rowPointers;
    std::vector<int> columnIndices;
    std::vector<double> values;
};

// Function to convert Matrix Market file to CSR format
CSRMatrix convertToCSR(const std::string& filename) {
    CSRMatrix csrMatrix;

    // Open the Matrix Market file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return csrMatrix;
    }

    // Variables to store matrix properties
    std::string line;
    int rows, cols, nnz;
    bool isPattern = false;

    // Read the Matrix Market header
    while (std::getline(file, line)) {
        if (line[0] != '%') {
            std::istringstream iss(line);
            iss >> rows >> cols >> nnz;
            if (line.find("pattern") != std::string::npos)
                isPattern = true;
            break;
        }
    }

    // Variables for CSR data structure
    csrMatrix.rowPointers.resize(rows + 1, 0);
    csrMatrix.columnIndices.resize(nnz, 0);
    if (!isPattern)
        csrMatrix.values.resize(nnz, 0.0);

    // Read the non-zero entries
    int row, col;
    double value = 1.0; // Default value if pattern-only
    int index = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> row >> col;
        if (!isPattern)
            iss >> value;

        row--; // Convert to 0-based indexing
        col--;

        csrMatrix.rowPointers[row + 1]++;
        csrMatrix.columnIndices[index] = col;
        if (!isPattern)
            csrMatrix.values[index] = value;

        index++;
    }

    // Calculate the row pointers
    for (int i = 2; i <= rows; i++) {
        csrMatrix.rowPointers[i] += csrMatrix.rowPointers[i - 1];
    }

    // Close the file
    file.close();

    return csrMatrix;
}

#endif  // CSR_CONVERSION_H
