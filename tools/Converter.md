# CSR Matrix Reader

This code provides functionality to read a matrix stored in the Matrix Market format and represent it using the Compressed Sparse Row (CSR) format. It also includes a function to print the CSR matrix representation.

## Dependencies

This code requires the following libraries:

- `iostream`: Provides input/output stream functionality.
- `fstream`: Provides file stream functionality for reading files.
- `vector`: Provides dynamic array functionality.

## Data Structure

The code defines a data structure called `CSRMatrix` to represent the matrix in CSR format. It consists of the following components:

- `rowPointers`: A vector of integers that stores the indices where each row starts in the `columnIndices` and `values` vectors.
- `columnIndices`: A vector of integers that stores the column indices of the non-zero elements.
- `values`: A vector of doubles that stores the values of the non-zero elements.

## Functions

### `bool readMatrixMarket(const std::string& filename, CSRMatrix& matrix)`

This function reads a matrix from a Matrix Market file and converts it to the CSR format.

**Parameters:**
- `filename`: A string containing the path to the Matrix Market file.
- `matrix`: A reference to a `CSRMatrix` object where the resulting CSR matrix will be stored.

**Returns:**
- `true` if the matrix is successfully read and converted.
- `false` if there was an error opening the file or reading the matrix.

### `void printCSRMatrix(const CSRMatrix& matrix)`

This function prints the contents of the CSR matrix.

**Parameters:**
- `matrix`: A constant reference to a `CSRMatrix` object to be printed.

### `int main(int argc, char* argv[])`

The main function of the program. It reads the command-line arguments, extracts the input filename, reads the matrix using the `readMatrixMarket` function, and prints the CSR matrix representation using the `printCSRMatrix` function.

**Command-line Arguments:**
- `-i <input_file.mtx>`: Specifies the input Matrix Market file.

## Usage

To use this code, compile it and run the resulting executable with the following command-line arguments:


### `./program -i <input_file.mtx>)`


- Replace `./program` with the name of the compiled executable.
- Replace `<input_file.mtx>` with the path to the input Matrix Market file.

The program will read the matrix from the file, convert it to CSR format, and print the CSR matrix representation to the console.

Note: Make sure to have the necessary Matrix Market file available for input.
