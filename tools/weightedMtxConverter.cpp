#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>
#include <getopt.h>

struct Edge {
    int source;
    int target;
};

void generateWeightedMTX(const std::string& inputFilename, const std::string& outputFilename, double minWeight, double maxWeight) {
    std::ifstream inputFile(inputFilename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }

    std::ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        std::cerr << "Error creating output file." << std::endl;
        inputFile.close();
        return;
    }

    std::string line;
    std::getline(inputFile, line); // Skip the header line

    int numRows, numCols, numEdges;
    inputFile >> numRows >> numCols >> numEdges;

    std::vector<Edge> edges;
    edges.reserve(numEdges);

    for (int i = 0; i < numEdges; ++i) {
        int row, col;
        inputFile >> row >> col;
        edges.push_back({row - 1, col - 1});
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<> dist(minWeight, maxWeight);

    //outputFile << "%%MatrixMarket matrix coordinate real general" << std::endl;
    //outputFile << "%%==================================================================================" << std::endl;
    outputFile << numRows << " " << numCols << " " << numEdges << std::endl;

    for (const auto& edge : edges) {
        double weight = dist(rng);
        outputFile << edge.source + 1 << " " << edge.target + 1 << " " << weight << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Weighted MTX file generated successfully." << std::endl;
}

int main(int argc, char* argv[]) {
    std::string inputFilename;
    std::string outputFilename;
    double minWeight = 0.1;
    double maxWeight = 1.0;

    int option;
    while ((option = getopt(argc, argv, "i:o:s:e:")) != -1) {
        switch (option) {
            case 'i':
                inputFilename = optarg;
                break;
            case 'o':
                outputFilename = optarg;
                break;
            case 's':
                minWeight = std::atof(optarg);
                break;
            case 'e':
                maxWeight = std::atof(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file> -o <output_file> -s <min_weight> -e <max_weight>" << std::endl;
                return 1;
        }
    }

    if (inputFilename.empty() || outputFilename.empty()) {
        std::cerr << "Input and output file names are required." << std::endl;
        return 1;
    }

    generateWeightedMTX(inputFilename, outputFilename, minWeight, maxWeight);

    return 0;
}
