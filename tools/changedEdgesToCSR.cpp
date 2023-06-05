#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// Struct to represent an edge in the modified CSR matrix
struct ModifiedEdge {
    int source;
    int destination;
    double weight;
    int changeType; // +1 for insertion, -1 for deletion
};

// Function to parse the output file and extract the modified edges
void parseOutputFile(const std::string& outputFile, std::vector<ModifiedEdge>& modifiedEdges) {
    std::ifstream file(outputFile);
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') {
            continue; // Skip empty lines and comment lines
        }

        std::istringstream iss(line);
        ModifiedEdge edge;
        iss >> edge.source >> edge.destination >> edge.weight >> edge.changeType;
        modifiedEdges.push_back(edge);
    }

    file.close();
}

// Function to print the modified CSR matrix
void printModifiedCSR(const std::vector<ModifiedEdge>& modifiedEdges) {
    std::cout << "Modified CSR Matrix:" << std::endl;
    for (const auto& edge : modifiedEdges) {
        std::cout << edge.source << "\t" << edge.destination << "\t" << edge.weight << "\t" << edge.changeType << std::endl;
    }
}

int main(int argc, char** argv) {
    // Check the command-line arguments
    if (argc != 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "-o") {
        std::cout << "Usage: ./modifiedCSRSpecifier -i <input_file> -o <output_file>\n";
        return 1;
    }

    // Parse the command-line arguments
    std::string inputFile = argv[2];
    std::string outputFile = argv[4];

    // Parse the output file and extract the modified edges
    std::vector<ModifiedEdge> modifiedEdges;
    parseOutputFile(outputFile, modifiedEdges);

    // Print the modified CSR matrix
    printModifiedCSR(modifiedEdges);

    return 0;
}
