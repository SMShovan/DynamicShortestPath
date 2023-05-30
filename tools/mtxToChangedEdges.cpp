#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cstring>

// Structure to represent an edge
struct Edge {
    int source;
    int destination;
    double weight;
};

// Function to parse Matrix Market file and populate the edges
std::vector<Edge> parseMatrixMarket(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Edge> edges;
    std::string line;

    // Skip comment lines
    while (getline(file, line) && line[0] == '%') {}

    int numNodes, numEdges;
    sscanf(line.c_str(), "%d %d", &numNodes, &numEdges);

    edges.reserve(numEdges);

    // Read edge data
    for (int i = 0; i < numEdges; i++) {
        Edge edge;
        getline(file, line);
        sscanf(line.c_str(), "%d %d %lf", &edge.source, &edge.destination, &edge.weight);
        edges.push_back(edge);
    }

    file.close();
    return edges;
}

// Function to write the modified edges to the output file
void writeModifiedEdges(const std::string& filename, const std::vector<Edge>& modifiedEdges) {
    std::ofstream file(filename);

    // Copy the comment lines from the input file
    std::ifstream inputFile("input.mtx");
    std::string line;
    while (getline(inputFile, line) && line[0] == '%') {
        file << line << "\n";
    }
    inputFile.close();

    file << modifiedEdges.size() << " " << 5 << " " << modifiedEdges.size() << "\n";

    for (const Edge& edge : modifiedEdges) {
        file << edge.source << " " << edge.destination << " " << edge.weight;

        if (edge.weight > 0) {
            file << " +1\n";
        } else {
            file << " -1\n";
        }
    }

    file.close();
}

// Function to randomly generate a positive edge weight
double generatePositiveWeight() {
    static std::mt19937 generator(time(0)); // Seed the generator with current time
    std::uniform_real_distribution<double> distribution(0.1, 10.0);
    return distribution(generator);
}




// Modify the edges based on the given insertion and deletion percentages
std::vector<Edge> modifyEdges(const std::vector<Edge>& edges, double insertionPercentage, double deletionPercentage, int totalChanges) {
    std::vector<Edge> modifiedEdges;

    int numInsertions = insertionPercentage * totalChanges;
    int numDeletions = deletionPercentage * totalChanges;

    // Delete modified edges
    std::vector<Edge> availableEdges = edges;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::shuffle(availableEdges.begin(), availableEdges.end(), generator);
    
    for (int i = 0; i < numDeletions && i < availableEdges.size(); i++) {
        Edge edge = availableEdges[i];
        edge.weight = -edge.weight;
        modifiedEdges.push_back(edge);
    }

    // Insert modified edges
    for (int i = 0; i < numInsertions; i++) {
        int source = rand() % edges.size() + 1;
        int destination = rand() % edges.size() + 1;
        double weight = generatePositiveWeight();
        modifiedEdges.push_back({source, destination, weight});
    }

    return modifiedEdges;
}



int main(int argc, char** argv) {
    // Parse command line arguments
    if (argc != 11) {
        std::cout << "Invalid number of arguments. Usage: ./script -f <input_file> -i <insertion_percentage> -d <deletion_percentage> -c <total_changes> -o <output_file>\n";
        return 1;
    }

    std::string inputFile;
    double insertionPercentage = 0.0;
    double deletionPercentage = 0.0;
    int totalChanges = 0;
    std::string outputFile;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0) {
            inputFile = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0) {
            insertionPercentage = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0) {
            deletionPercentage = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "-c") == 0) {
            totalChanges = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0) {
            outputFile = argv[++i];
        }
    }

    // Parse the input Matrix Market file
    std::vector<Edge> edges = parseMatrixMarket(inputFile);

    // Modify the edges based on the given insertion and deletion percentages
    std::vector<Edge> modifiedEdges = modifyEdges(edges, insertionPercentage, deletionPercentage, totalChanges);

    // Write the modified edges to the output file
    writeModifiedEdges(outputFile, modifiedEdges);

    std::cout << "Modified edges written to " << outputFile << "\n";

    return 0;
}
