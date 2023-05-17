#include <iostream>
#include <fstream>
#include <vector>

struct CSRMatrix {
    std::vector<int> rowPointers;
    std::vector<int> columnIndices;
    std::vector<double> values;
};

bool readMatrixMarket(const std::string& filename, CSRMatrix& matrix) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;
    int numRows, numCols, numNonZeros;

    // Skip comment lines
    while (std::getline(file, line) && line[0] == '%');

    // Read matrix metadata
    sscanf(line.c_str(), "%d %d %d", &numRows, &numCols, &numNonZeros);

    matrix.rowPointers.resize(numRows + 1, 0);

    // Read non-zero entries
    int row, col;
    double value;
    for (int i = 0; i < numNonZeros; i++) {
        file >> row >> col >> value;
        matrix.values.push_back(value);
        matrix.columnIndices.push_back(col - 1);
        matrix.rowPointers[row]++;
    }

    // Calculate row pointers
    int prefixSum = 0;
    for (int i = 0; i <= numRows; i++) {
        int count = matrix.rowPointers[i];
        matrix.rowPointers[i] = prefixSum;
        prefixSum += count;
    }

    file.close();
    return true;
}

void printCSRMatrix(const CSRMatrix& matrix) {
    std::cout << "Row Pointers: ";
    for (int i = 0; i < matrix.rowPointers.size(); i++) {
        std::cout << matrix.rowPointers[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Column Indices: ";
    for (int i = 0; i < matrix.columnIndices.size(); i++) {
        std::cout << matrix.columnIndices[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Values: ";
    for (int i = 0; i < matrix.values.size(); i++) {
        std::cout << matrix.values[i] << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: ./program -i <input_file.mtx>" << std::endl;
        return 1;
    }

    std::string inputFilename;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-i") {
            if (i + 1 < argc) {
                inputFilename = argv[i + 1];
                break;
            }
        }
    }

    if (inputFilename.empty()) {
        std::cout << "Input filename not provided." << std::endl;
        return 1;
    }

    CSRMatrix matrix;
    if (!readMatrixMarket(inputFilename, matrix)) {
        return 1;
    }

    std::cout << "CSR Matrix Representation:" << std::endl;
    printCSRMatrix(matrix);

    return 0;
}
