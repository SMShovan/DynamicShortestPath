# mtxToChangedEdges.cpp

The "mtxToChangedEdges.cpp" script takes an input Matrix Market (.mtx) file and generates an output file that contains the newly inserted and deleted edges based on the specified percentage of insertion and deletion. The output file is also in the Matrix Market format.

## Usage

The script is executed with the following command-line arguments:

```plaintext 
%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This is a comment line
%
%=================================================================================
4 4 6
1 1 2.0
1 2 3.5
2 3 4.2
3 4 1.8
3 3 5.1
```


We want to generate an output file called "output.mtx" with 30% insertions, 20% deletions, and a total of 10 changes. To achieve this, we run the script with the following command:

```
./mtxToChangedEdges -f input.mtx -i 0.3 -d 0.2 -c 10 -o output.mtx
```


The script will process the input file and generate the output file "output.mtx" with the modified edges. The output file will have the following structure:

```plaintext

%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This is a comment line
%
%=================================================================================
4 5 10
1 1 2.0 +1
1 2 3.5 +1
2 3 4.2 +1
3 4 1.8 +1
3 3 5.1 +1
4 2 2.7 +1
1 1 2.0 -1
3 3 5.1 -1
3 4 1.8 -1
4 2 2.7 -1
```


The output file retains the same comment lines from the input file. It has 4 rows, 5 columns, and a total of 10 modified edges (including both insertions and deletions). The modified edges are denoted with a `+1` or `-1` in the rightmost column, indicating insertion or deletion, respectively. The edge weights for inserted edges are randomly assigned positive values, while the edge weights for deleted edges preserve their original values.

Please note that the specific edges modified and the randomly generated values may differ based on the execution of the script.

---

This completes the documentation for the "mtxToChangedEdges.cpp" script. It allows you to generate an output file with the newly inserted and deleted edges based on the specified percentage of insertion and deletion. The script follows a step-by-step process to parse the input file, modify the edges, and write the modified edges to the output file. Feel free to copy this documentation and use it as needed.

Let me know if you have any further questions!
