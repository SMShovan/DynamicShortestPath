#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <queue>
#include <cstdlib>

using namespace std;

// Function to check if a graph is connected using BFS.
bool isConnected(const vector<vector<int>>& adjacency_matrix, int vertices) {
    vector<bool> visited(vertices + 1, false);
    queue<int> q;
    int count = 0;

    q.push(1);
    visited[1] = true;
    count++;

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        for (int i = 1; i <= vertices; i++) {
            if (adjacency_matrix[current][i] != 0 && !visited[i]) {
                q.push(i);
                visited[i] = true;
                count++;
            }
        }
    }

    return count == vertices;
}

int main(int argc, char* argv[]) {
    if (argc != 11) {
        cout << "Usage: ./random_graph -v <vertices> -e <edges> -b <min_range> -t <max_range> -o <output_file>\n";
        return 1;
    }

    int vertices, edges, min_range, max_range;
    string output_file;

    for (int i = 1; i < argc; i += 2) {
        string option = argv[i];
        string value = argv[i + 1];

        if (option == "-v")
            vertices = atoi(value.c_str());
        else if (option == "-e")
            edges = atoi(value.c_str());
        else if (option == "-b")
            min_range = atoi(value.c_str());
        else if (option == "-t")
            max_range = atoi(value.c_str());
        else if (option == "-o")
            output_file = value;
        else {
            cout << "Invalid option: " << option << "\n";
            return 1;
        }
    }

    vector<vector<int>> adjacency_matrix(vertices + 1, vector<int>(vertices + 1, 0));

    // Generate a random graph with the given number of vertices and edges.
    default_random_engine generator;
    uniform_int_distribution<int> distribution(min_range, max_range);

    for (int i = 0; i < edges; i++) {
        int u = rand() % vertices + 1;
        int v = rand() % vertices + 1;

        // Make sure that the graph is connected by adding an edge between u and v if they are not already connected.
        if (u != v && adjacency_matrix[u][v] == 0) {
            adjacency_matrix[u][v] = distribution(generator);
            adjacency_matrix[v][u] = distribution(generator);
        }
    }

    // Check if the graph is connected, if not, add additional edges to ensure connectivity.
    while (!isConnected(adjacency_matrix, vertices)) {
        int u = rand() % vertices + 1;
        int v = rand() % vertices + 1;

        if (u != v && adjacency_matrix[u][v] == 0) {
            adjacency_matrix[u][v] = distribution(generator);
            adjacency_matrix[v][u] = distribution(generator);
        }
    }

    // Write the adjacency matrix to the specified output file.
    ofstream outfile(output_file);
    outfile << vertices << " " << vertices << " " << edges << "\n";

    // Iterate over the rows of the adjacency matrix and add the nonzero elements to the MTX file.
    for (int i = 1; i <= vertices; i++) {
        for (int j = 1; j <= vertices; j++) {
            if (adjacency_matrix[i][j] != 0) {
                outfile << i << " " << j << " " << adjacency_matrix[i][j] << "\n";
            }
        }
    }

    outfile.close();

    cout << "Random graph generated and saved to '" << output_file << "'." << endl;

    return 0;
}

//clang++ -std=c++11 randomConnectedGraph.cpp -o randomConnectedGraph && ./randomConnectedGraph -v 5 -e 10 -b 1 -t 10 -o randomgraph.mtx
