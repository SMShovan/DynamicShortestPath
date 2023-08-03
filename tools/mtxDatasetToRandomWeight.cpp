#include <iostream>
#include <fstream>
#include <vector>
#include <random>

struct Edge {
    int row;
    int col;
    int weight;
};
int num_nodes;

std::vector<Edge> read_mtx_file(const std::string& file_path) {
    std::vector<Edge> edges;
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open the file " << file_path << std::endl;
        return edges;
    }

    std::string line;
    int i =0; 
    while (std::getline(file, line)) {
        if (line[0] == '%')
            continue;
        
        int vertex1,vertex2, edgeCount;
        if (sscanf(line.c_str(), "%d %d %d", &vertex1, &vertex2, &edgeCount) == 3) {
            num_nodes = vertex1;
        }
        

        int row, col;
        if (sscanf(line.c_str(), "%d %d", &row, &col) == 2) {
            i++; 
            if (i == 1)
                continue;
            edges.push_back({row, col, 0});
        }
    }

    file.close();
    return edges;
}

std::vector<Edge> convert_to_weighted_graph(const std::vector<Edge>& edges, int start_range, int end_range) {
    std::vector<Edge> weighted_edges;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(start_range, end_range);

    for (const auto& edge : edges) {
        int weight = dist(gen);
        weighted_edges.push_back({edge.row, edge.col, weight});
    }

    return weighted_edges;
}



void save_weighted_graph_to_mtx(const std::vector<Edge>& weighted_edges, const std::string& output_file_path) {
    std::ofstream file(output_file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to create the output file " << output_file_path << std::endl;
        return;
    }

    //file << "%%MatrixMarket matrix coordinate real symmetric\n";
    //file << "% Weighted graph with " << num_nodes << " nodes and " << weighted_edges.size() << " edges\n";
    file << num_nodes << " " << num_nodes << " " << weighted_edges.size() << "\n";

    for (const auto& edge : weighted_edges) {
        file << edge.row << " " << edge.col << " " << edge.weight << "\n";
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -s <start_range> -e <end_range>" << std::endl;
        return 1;
    }

    std::string input_file_path;
    std::string output_file_path;
    int start_range = 1;
    int end_range = 10;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-i") {
            if (i + 1 < argc) {
                input_file_path = argv[i + 1];
            } else {
                std::cerr << "Error: Input file path not specified." << std::endl;
                return 1;
            }
        } else if (arg == "-o") {
            if (i + 1 < argc) {
                output_file_path = argv[i + 1];
            } else {
                std::cerr << "Error: Output file path not specified." << std::endl;
                return 1;
            }
        } else if (arg == "-s") {
            if (i + 1 < argc) {
                start_range = std::atoi(argv[i + 1]);
                if (start_range < 1) {
                    std::cerr << "Error: Starting range must be a positive integer." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: Starting range not specified." << std::endl;
                return 1;
            }
        } else if (arg == "-e") {
            if (i + 1 < argc) {
                end_range = std::atoi(argv[i + 1]);
                if (end_range < start_range) {
                    std::cerr << "Error: Ending range must be greater than or equal to the starting range." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: Ending range not specified." << std::endl;
                return 1;
            }
        }
    }

    if (input_file_path.empty() || output_file_path.empty()) {
        std::cerr << "Error: Both input and output file paths must be specified." << std::endl;
        return 1;
    }

    std::vector<Edge> edges = read_mtx_file(input_file_path);
    std::vector<Edge> weighted_edges = convert_to_weighted_graph(edges, start_range, end_range);
    save_weighted_graph_to_mtx(weighted_edges, output_file_path);

    return 0;
}
