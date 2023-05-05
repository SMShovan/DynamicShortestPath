#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

struct MCSR {
    std::vector<int> row_ptr;    // row pointer array
    std::vector<int> block_ptr;  // block pointer array
    std::vector<int> block_size; // block size array
    std::vector<int> col_index;  // column index array
    std::vector<int> data;       // data array
};

void print_mcsr(const MCSR& mcsr_data) {
    std::cout << "Modified CSR format:\n";
    std::cout << "Num rows: " << mcsr_data.row_ptr.size() - 1 << "\n";
    std::cout << "Num nonzeros: " << mcsr_data.col_index.size() << "\n";
    std::cout << "Row pointer array:\n";
    for (int i = 0; i < mcsr_data.row_ptr.size(); i++) {
        std::cout << mcsr_data.row_ptr[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Block pointer array:\n";
    for (int i = 0; i < mcsr_data.block_ptr.size(); i++) {
        std::cout << mcsr_data.block_ptr[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Block size array:\n";
    for (int i = 0; i < mcsr_data.block_size.size(); i++) {
        std::cout << mcsr_data.block_size[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Column index array:\n";
    for (int i = 0; i < mcsr_data.col_index.size(); i++) {
        std::cout << mcsr_data.col_index[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Data array:\n";
    for (int i = 0; i < mcsr_data.data.size(); i++) {
        std::cout << mcsr_data.data[i] << " ";
    }
    std::cout << "\n";
}

MCSR mtx_to_mcsr(std::string mtx_file) {
    // Open input file
    std::ifstream file(mtx_file);

    // Read header information
    std::string line;
    int num_rows, num_cols, num_nonzeros;
    getline(file, line);
    getline(file, line);
    getline(file, line);
    std::istringstream iss(line);
    iss >> num_rows >> num_cols >> num_nonzeros;

    // Initialize MCSR data structure
    MCSR mcsr;
    mcsr.row_ptr.resize(num_rows + 1, 0); // extra entry for total number of nonzeros
    mcsr.block_ptr.push_back(0);          // first block starts at index 0
    mcsr.col_index.reserve(num_nonzeros);
    mcsr.data.reserve(num_nonzeros);

    // Read non-zero entries from file
    int row, col;
    while (getline(file, line)) {
        std::istringstream iss(line);
        iss >> row >> col;
        // Convert 1-indexed to 0-indexed
        row--;
        col--;
        mcsr.col_index.push_back(col);
        mcsr.data.push_back(1); // assuming pattern matrix
        // Update row pointer array
        mcsr.row_ptr[row + 1]++;
    }

        // Compute row pointer array
    for (int i = 1; i <= num_rows; i++) {
        mcsr.row_ptr[i] += mcsr.row_ptr[i - 1];
    }
    mcsr.row_ptr.push_back(num_nonzeros);

    // Compute block pointer and size arrays
    int block_size = num_rows / 2 + 1;
    int num_blocks = num_rows / block_size + (num_rows % block_size == 0 ? 0 : 1);
    mcsr.block_size.resize(num_blocks);
    for (int i = 0; i < num_blocks - 1; i++) {
        mcsr.block_size[i] = block_size;
        mcsr.block_ptr.push_back(mcsr.block_ptr[i] + block_size * (mcsr.row_ptr[(i+1)*block_size] - mcsr.row_ptr[i*block_size]));
    }
    mcsr.block_size[num_blocks - 1] = num_rows - (num_blocks - 1) * block_size;
    mcsr.block_ptr.push_back(num_nonzeros);

    return mcsr;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_file.mtx\n";
        return 1;
    }

    std::string input_file = argv[1];
    MCSR mcsr_data = mtx_to_mcsr(input_file);

    print_mcsr(mcsr_data);

    return 0;
}
