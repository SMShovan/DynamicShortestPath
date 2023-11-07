#include <bits/stdc++.h>
struct Edge {
    int source;
    int destination;
    double weight;

    Edge(int src, int dest, double w) : source(src), destination(dest), weight(w) {}
};
using Graph = std::unordered_map<int, std::unordered_map<int, int>>;

int& access_with_default(Graph& g, int key1, int key2, int default_value = 4) {
    if(g[key1].find(key2) == g[key1].end()) {
        g[key1][key2] = default_value;
    }
    return g[key1][key2];
}

class Tree {
public:
    virtual const std::vector<int>& getParents() const = 0;  // Pure virtual method
    virtual ~Tree() = default;  // Virtual destructor
};

class Tree1 : public Tree {
    std::vector<int> parents;
public:
    Tree1(const std::vector<int>& p) : parents(p) {}
    const std::vector<int>& getParents() const override {
        return parents;
    }
};

class Tree2 : public Tree {
    std::vector<int> parents;
public:
    Tree2(const std::vector<int>& p) : parents(p) {}
    const std::vector<int>& getParents() const override {
        return parents;
    }
};

class Tree3 : public Tree {
    std::vector<int> parents;
public:
    Tree3(const std::vector<int>& p) : parents(p) {}
    const std::vector<int>& getParents() const override {
        return parents;
    }
};

void addEdge(Graph& graph, int u, int v, double pref) {
    access_with_default(graph, u, v) = access_with_default(graph, u, v) - 1.0/pref;
}

Graph constructGraph(const std::vector<Tree*>& trees, const std::vector<double>& Pref) {
    Graph graph;

    for (size_t index = 0; index < trees.size(); ++index) {
        const auto& tree = trees[index];
        const std::vector<int>& parents = tree->getParents();

        for (int i = 1; i < parents.size(); ++i) {
            if (parents[i] != 0) {
                addEdge(graph, parents[i], i, Pref[index]);
            }
        }
    }
    return graph;
}

bool bellmanFord(const Graph& graph, int source, std::unordered_map<int, int>& distances, std::vector<int>& newParent) {
    int V = graph.size();
    newParent.assign(V + 1, -1);  // Initialize with -1 (no predecessor)

    // Initialize distances
    for (const auto& [node, _] : graph) {
        distances[node+1] = std::numeric_limits<int>::max();
    }
    distances[source] = 0;

    // Relax edges |V| - 1 times
    for (size_t i = 1; i <= V; ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            auto it = std::next(graph.begin(), j);
            int u = it->first;
            const auto& neighbors = it->second;
            for (const auto& [v, weight] : neighbors) {
                if (distances[u] != std::numeric_limits<int>::max() && distances[u] + weight < distances[v]) {
   
                    if (distances[u] + weight < distances[v]) {
                        distances[v] = distances[u] + weight;
                        newParent[v] = u;
                    }

                }
            }
        }
    }

    // Check for negative-weight cycles
    for (const auto& [u, neighbors] : graph) {
        for (const auto& [v, weight] : neighbors) {
            if (distances[u] != std::numeric_limits<int>::max() && distances[u] + weight < distances[v]) {
                return false;  
            }
        }
    }

    return true;  
}

bool operator==(const Edge& lhs, const Edge& rhs) {
    return lhs.source == rhs.source && lhs.destination == rhs.destination && lhs.weight == rhs.weight;
}

// The following function reads the graph or changed edges files, however, the DT name should be changed as it is not storing in DT format. We only calculate predecessor if the file contains the graph.
std::vector<std::vector<Edge>> convertToDT(std::ifstream& inputFile, bool isGraph, std::vector<std::vector<Edge>>& predecessor) {
    std::string line;
    int numRows, numCols, numNonZero;
    if (!std::getline(inputFile, line)) {
        std::cerr << "Error: Unable to read the header line from the input file." << std::endl;
        exit(1);
    }
    std::istringstream headerStream(line);
    if (!(headerStream >> numRows >> numCols >> numNonZero)) {
        std::cerr << "Error: Invalid header format in the input file." << std::endl;
        exit(1);
    }
    std::vector<std::vector<Edge>> DTMatrix(numRows);
    int lineCount = 1;
    while (std::getline(inputFile, line)) {
        lineCount++;
        std::istringstream lineStream(line);
        int row, col;
        double value;
        if (!(lineStream >> row >> col >> value)) {
            std::cerr << "Error: Invalid line format at line " << lineCount << " in the input file." << std::endl;
            exit(1);
        }
        if (row < 1 || row > numRows || col < 1 || col > numCols) {
            std::cerr << "Error: Invalid vertex indices at line " << lineCount << " in the input file." << std::endl;
            exit(1);
        }
        DTMatrix[row - 1].emplace_back(row - 1, col - 1, value);
        if(isGraph)
        {
            predecessor.resize(numCols);
            predecessor[col - 1].emplace_back(row, col, value);
        }
    }    
    return DTMatrix;
}

std::vector<std::vector<int>> dijkstra(const std::vector<std::vector<Edge>>& graphDT, int sourceNode, std::vector<double>& shortestDist, std::vector<int>& parentList) {
    int numNodes = graphDT.size();

    std::vector<bool> visited(numNodes, false);

    shortestDist[sourceNode - 1] = 0;

    for (int i = 0; i < numNodes - 1; ++i) {
        int minDistNode = -1;
        double minDist = std::numeric_limits<double>::infinity();
        for (int j = 0; j < numNodes; ++j) {
            if (!visited[j] && shortestDist[j] < minDist) {
                minDist = shortestDist[j];
                minDistNode = j;
            }
        }
        if (minDistNode == -1) {
            break;
        }
        visited[minDistNode] = true;

        for (const Edge& edge : graphDT[minDistNode]) {
            int neighbor = edge.destination;
            double weight = edge.weight;
            if (!visited[neighbor] && shortestDist[minDistNode] + weight < shortestDist[neighbor]) {
                shortestDist[neighbor] = shortestDist[minDistNode] + weight;
            }
        }
    }
    
    // Build the shortest path tree based on the shortest distances
    std::vector<std::vector<int>> ssspTree(numNodes);
    std::vector<bool> cycleCheck (numNodes, false);
    for (int i = 0; i < numNodes; ++i) {
        if (shortestDist[i] != std::numeric_limits<double>::infinity()) {
            int parent = i + 1;
            for (const Edge& edge : graphDT[i]) {
                int child = edge.destination + 1;
                if (shortestDist[child - 1] == shortestDist[i] + edge.weight && !cycleCheck[child - 1]) {
                    ssspTree[parent - 1].push_back(child);
                    cycleCheck[child - 1 ] = true;
                    parentList[child] = parent;
                }
            }
        }
    }
    return ssspTree;
}

void printShortestPathTree(const std::vector<std::pair<int, std::vector<int>>>& parentChildSSP) {
    std::cout << "Shortest Path Tree:\n";
    for (const auto& node : parentChildSSP) {
        std::cout << "Node " << node.first << ": ";
        for (int child : node.second) {
            std::cout << child << " ";
        }
        std::cout << "\n";
    }
}

void markSubtreeAffected(std::vector<std::pair<int, std::vector<int>>>& parentChildSSP, std::vector<double>& shortestDist, std::vector<bool>& affectedNodes, std::queue<int>& affectedNodesQueue, std::vector<bool>& affectedDelNodes, int node) {

    for (int child : parentChildSSP[node].second) {
        if( !affectedDelNodes[child - 1])
            affectedDelNodes[child - 1] = true;
    }
}

std::vector<std::pair<int, std::vector<int>>> convertToParentChildSSP(const std::vector<std::vector<int>>& ssspTree) {
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP(ssspTree.size());

    for (int i = 0; i < ssspTree.size(); ++i) {
        parentChildSSP[i].first = i + 1;  
        parentChildSSP[i].second = ssspTree[i];  
    }
    return parentChildSSP;
}

/*
1. ssspTree is a subset (0-indexed) of graphDT having the edges belongs to shortest path. However, ssspTree contains a set of pairs with parent (1-indexed) and children (1-indexed).
2. graphDT is the whole graph containg a set of vector of edges. It is 0-indexed. i index contains all the outgoing edges (also 0-indexed) from vertex (i+1). The graph is assumed to have vertex id > 0. 
3. shortestDist is the set of all vertices (0-indexed) from the source node assumed to be 1 (in 1-index).
4. parentList is the set of all parents (1-indexed) given the child node (1-indexed) of the ssspTree.
5. Global variable Predecessor is similar to graphDT, however, instead of storing into rows, it stores in column index to find all incident vertices. Again the column index is 0-indexed, but edges are 1-indexed.
*/

void updateShortestPath(std::vector<std::pair<int, std::vector<int>>>& ssspTree, std::vector<std::vector<Edge>>& graphDT, const std::vector<std::vector<Edge>>& changedEdgesDT, std::vector<double>& shortestDist, std::vector<int>& parentList, std::vector<std::vector<Edge>>& predecessor) {
    std::vector<Edge> insertedEdges;
    std::vector<Edge> deletedEdges;
    std::vector<bool> affectedNodes(ssspTree.size(), false);
    std::vector<bool> affectedDelNodes(ssspTree.size(), false);
    std::queue<int> affectedNodesQueue;

    // Identify inserted and deleted edges
    for (const std::vector<Edge>& row : changedEdgesDT) {
        for (const Edge& edge : row) {
            if (edge.weight > 0) {
                insertedEdges.push_back(edge);
            } else if (edge.weight < 0) {
                deletedEdges.push_back(edge);
            }
        }
    }

    for (const Edge& edge : insertedEdges) {
        
        int x,y; 
        x = edge.source;
        y = edge.destination;

        graphDT[x].push_back(Edge(x, y, edge.weight));
        predecessor[y].emplace_back(x, y, edge.weight);
        if (shortestDist[y] > shortestDist[x] + edge.weight) {
            shortestDist[y] = shortestDist[x] + edge.weight;

            int oldParent = parentList[y + 1];

            for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
            {
                if (ssspTree[oldParent - 1].second[j] == y + 1 )
                {
                    ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                    break;
                }
            }
            parentList[y+1] = x + 1;
            ssspTree[x].second.push_back(y+1);
            affectedNodes[y] = true; 
            affectedNodesQueue.push(y);
        }
    }

    for (const Edge& edge : deletedEdges) {
        for (int i = 0; i < graphDT[edge.source].size(); ++i) {
        if (graphDT[edge.source][i].source == edge.source && graphDT[edge.source][i].destination == edge.destination) {
            graphDT[edge.source].erase(graphDT[edge.source].begin() + i);
            break;
        }
        }
        bool inTree = false;
        for (auto& preEdge : predecessor[edge.destination]) 
        if (edge.source + 1 == preEdge.source) {
            
            // Find the iterator pointing to the element
            auto it = std::find(predecessor[edge.destination].begin(), predecessor[edge.destination].end(), preEdge);
            
            for (int i = 0; i < ssspTree[edge.source].second.size() ; i++)
            {
                if (ssspTree[edge.source].second[i] == edge.destination + 1)
                inTree = true;
            }

            if (it != predecessor[edge.destination].end()) {
                predecessor[edge.destination].erase(it);
            }
        }


        if (!inTree)
            continue;

        affectedNodes[edge.destination] = true;
        affectedNodesQueue.push(edge.destination);
        affectedDelNodes[edge.destination] = true;
        shortestDist[edge.destination] = std::numeric_limits<double>::infinity();

        markSubtreeAffected(ssspTree, shortestDist, affectedNodes, affectedNodesQueue, affectedDelNodes, edge.destination);

        int n = edge.destination; 
        int newDistance = std::numeric_limits<double>::infinity();
        int newParentIndex = -1; 

        for (int i = 0; i < predecessor[n].size(); i++)
        {
            if (shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight < newDistance)
            {
                newDistance = shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight;

                newParentIndex = predecessor[n][i].source - 1;
            } 
        }
        int oldParent = parentList[edge.destination + 1];
        
        if (newParentIndex == -1)
        {
            parentList[n + 1] = -1; 
            shortestDist[n] = std::numeric_limits<double>::infinity();
        }   
        else 
        {    
            shortestDist[n] = newDistance; 
            ssspTree[newParentIndex].second.push_back(n+1);   
            parentList[edge.destination + 1] = newParentIndex + 1;   
        }
            for (int j = 0; j < ssspTree[oldParent - 1].second.size() && oldParent >= 1 ; j++ )
            {
                
                if (ssspTree[oldParent - 1].second[j] == edge.destination + 1 )
                {
                    ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                    break;
                }
            }
    }
        
    while(!affectedNodesQueue.empty()){
        int v = affectedNodesQueue.front();
        affectedNodesQueue.pop();
        affectedNodes[v] = false;  

        for (const Edge& edge : graphDT[v] ) {
            int n = edge.destination;

            // for deletion
            if ( shortestDist[v] == std::numeric_limits<double>::infinity() ){
                if ( n + 1 == 1)
                    continue;

                int newDistance = std::numeric_limits<double>::infinity();
                int newParentIndex = -1; 
                for (int i = 0; i < predecessor[n].size(); i++)
                {
                    if (shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight < newDistance)
                    {
                        newDistance = shortestDist[predecessor[n][i].source - 1] + predecessor[n][i].weight;
                        newParentIndex = predecessor[n][i].source - 1;
                    } 
                }

                int oldParent = parentList[n + 1];

                if (newParentIndex == -1)
                {
                    parentList[n + 1] = -1; 
                    shortestDist[n] = std::numeric_limits<double>::infinity();
                }
                else 
                {                        
                    shortestDist[n] = newDistance; 
                    ssspTree[newParentIndex].second.push_back(n+1);  
                    parentList[n + 1] = newParentIndex + 1;     
                }    

                for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
                {
                    
                    if (ssspTree[oldParent - 1].second[j] == n + 1 )
                    {
                        ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                        break;
                    }
                }
                 
                affectedNodes[n] = true;  
                affectedNodesQueue.push(n);
            }

            // for insertion
            if (shortestDist[n] > shortestDist[v] + edge.weight && shortestDist[v] != std::numeric_limits<double>::infinity()) 
            {    
                shortestDist[n] = shortestDist[v] + edge.weight;
                int oldParent = parentList[n + 1];
                for (int j = 0; j < ssspTree[oldParent - 1].second.size(); j++ )
                {
                    if (ssspTree[oldParent - 1].second[j] == n + 1 )
                    {
                        ssspTree[oldParent - 1].second.erase(ssspTree[oldParent - 1].second.begin() + j);
                        break;
                    }
                }
                parentList[n + 1] = v + 1; 
                ssspTree[v].second.push_back(n+1); 
                affectedNodes[n] = true; 
                affectedNodesQueue.push(n);
            }
        }
    }
}

int main(int argc, char** argv) 
{      
    if (argc < 4) {
        std::cerr << "Usage: ./program -g <graph_file> -c <changed_edges_file> -s <source_node>\n";
        return 1;
    }
    std::string graphFile;
    std::string changedEdgesFile;
    int sourceNode;
    
    for (int i = 1; i < argc; i += 2) {
        std::string option(argv[i]);
        std::string argument(argv[i + 1]);

        if (option == "-g") {
            graphFile = argument;
        } else if (option == "-c") {
            changedEdgesFile = argument;
        } else if (option == "-s") {
            sourceNode = std::stoi(argument);
        } else {
            std::cerr << "Invalid option: " << option << "\n";
            return 1;
        }
    }
    std::ifstream graphInputFile(graphFile);
    if (!graphInputFile.is_open()) {
        std::cerr << "Unable to open graph file: " << graphFile << std::endl;
        return 1;
    }
    std::vector<std::vector<Edge>> predecessor;
    std::vector<std::vector<Edge>> graphDT = convertToDT(graphInputFile, true, predecessor);
    graphInputFile.close();
    std::vector<double> shortestDist(graphDT.size(), std::numeric_limits<double>::infinity());
    shortestDist[sourceNode - 1] = 0;
    std::vector<int> parent(graphDT.size() + 1, -1);
    std::vector<std::vector<int>> ssspTree = dijkstra(graphDT, sourceNode, shortestDist, parent);
    std::vector<std::pair<int, std::vector<int>>> parentChildSSP = convertToParentChildSSP(ssspTree);
    std::ifstream changedEdgesInputFile(changedEdgesFile);
    std::vector<std::vector<Edge>> changedEdgesDT = convertToDT(changedEdgesInputFile, false, predecessor);
    changedEdgesInputFile.close();

    updateShortestPath(parentChildSSP, graphDT, changedEdgesDT, shortestDist, parent, predecessor);

    printShortestPathTree(parentChildSSP);

    Tree1 t1(parent);
    Tree2 t2(parent);
    Tree3 t3(parent);


    std::vector<Tree*> trees = {&t1, &t2, &t3};
    std::vector<double> Pref = {1.0,1.0,1.0};
    Graph graph = constructGraph(trees, Pref);

    // Print the graph
    for (const auto& [node, neighbors] : graph) {
        std::cout << "Node " << node << " has directed edges to:\n";
        for (const auto& [neighbor, weight] : neighbors) {
            std::cout << "  Node " << neighbor << " with weight " << weight << "\n";
        }
    }
    std::unordered_map<int, int> distances;
    std::vector<int> newParent;
    int source = 1;  // Example source node
    bool success = bellmanFord(graph, source, distances, newParent);

    if (success) {
        std::cout << "Shortest distances from node " << source << ":\n";
        for (const auto& [node, distance] : distances) {
            std::cout << "Node " << node << ": " << distance << ", Parent: " << newParent[node] << "\n";
        }
    } else {
        std::cout << "The graph contains a negative-weight cycle.\n";
    }

    return 0;
}


// To Run: clang++ -std=c++11 SOSP.cpp -o program && ./program -g input_graph.mtx -c changed_edges.mtx -s 1      
