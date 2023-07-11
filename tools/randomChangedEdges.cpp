#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cstring>

using namespace std;

struct Edge {
    int u;
    int v;
    int weight;
};

bool isEdgeExists(const vector<Edge>& edges, int row, int col) {
    for (const Edge& edge : edges) {
        if ((edge.u == row && edge.v == col) || (edge.u == col && edge.v == row)) {
            return true;
        }
    }
    return false;
}

int main(int argc, char* argv[]) {
    if (argc != 9) {
        cout << "Usage: ./changeEdges -i <num_inserted_edges> -d <num_deleted_edges> -g <original_graph.mtx> -o <output_file>\n";
        return 1;
    }

    int numInsertedEdges = 0;
    int numDeletedEdges = 0;
    string originalGraphFile;
    string output_file;

    // Parse command-line arguments
    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0) {
            numInsertedEdges = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-d") == 0) {
            numDeletedEdges = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-g") == 0) {
            originalGraphFile = argv[i + 1];
        }
        else if (strcmp(argv[i], "-o") == 0) {
            output_file = argv[i + 1];
        }
        else {
            cout << "Invalid option: " << argv[i] << "\n";
            return 1;
        }
    }

    ifstream original_graph(originalGraphFile);
    if (!original_graph) {
        cout << "Error: Failed to open the original graph file.\n";
        return 1;
    }

    int numVertices, numEdges;
    original_graph >> numVertices >> numVertices >> numEdges;

    vector<Edge> edges;
    for (int i = 0; i < numEdges; i++) {
        Edge edge;
        original_graph >> edge.u >> edge.v >> edge.weight;
        edges.push_back(edge);
    }

    original_graph.close();

    ofstream outfile(output_file);
    outfile << numVertices << " " << numVertices << " " << numInsertedEdges + numDeletedEdges << "\n";

    default_random_engine generator;
    uniform_int_distribution<int> distributionWeight(1, 10);
    uniform_int_distribution<int> distributionNode(1, numVertices);

    // Generate inserted edges with random positive weights.
    for (int i = 0; i < numInsertedEdges; i++) {
        int row = distributionNode(generator);
        int col = distributionNode(generator);
        int weight = distributionWeight(generator);

        outfile << row << " " << col << " " << weight << "\n";
    }

    // Generate deleted edges from the original graph.
    int deletedEdgesCount = 0;
    while (deletedEdgesCount < numDeletedEdges) {
        int randomEdgeIndex = distributionNode(generator) % numEdges;
        const Edge& edge = edges[randomEdgeIndex];

        // Check if the randomly selected edge exists and has not been deleted already.
        if (isEdgeExists(edges, edge.u, edge.v) && edge.weight > 0) {
            outfile << edge.u << " " << edge.v << " -1\n";
            deletedEdgesCount++;
        }
    }

    outfile.close();

    cout << "Change edges generated and saved to '" << output_file << "'." << endl;

    return 0;
}

//clang++ -std=c++11 randomChangedEdges.cpp -o randomChangedEdges && ./randomChangedEdges -i 2 -d 2 -g randomgraph.mtx -o changeEdges.mtx