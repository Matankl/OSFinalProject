#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

using namespace std;

class WeightedGraph {
    int n;  // Number of vertices
    vector<vector<pair<int, int>>> adj;  // Adjacency list with weight (vertex, weight)

public:
    WeightedGraph(int vertices) : n(vertices) {
        adj.resize(n + 1);  // Resize adjacency list for 1-based indexing
    }

    // Function to add an edge with weight
    void addEdge(int u, int v, int weight) {
        adj[u].emplace_back(v, weight);  // Add edge u -> v with weight
        adj[v].emplace_back(u, weight);  // For undirected graphs
    }

    // Get the adjacency list of the graph
    const vector<pair<int, int>>& getAdjList(int v) const {
        return adj[v];
    }

    // Get the number of vertices
    int getNumVertices() const {
        return n;
    }
};
