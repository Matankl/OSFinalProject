#ifndef WEIGHTED_GRAPH_HPP
#define WEIGHTED_GRAPH_HPP

#include <vector>
#include <utility>

class WeightedGraph {
    int n;  // Number of vertices
    std::vector<std::vector<std::pair<int, int>>> adj;  // Adjacency list with weight (vertex, weight)

public:
    // Constructor to initialize the graph with the number of vertices
    WeightedGraph(int vertices);

    // Function to add an edge with weight
    void addEdge(int u, int v, int weight);

    // Get the adjacency list of a vertex
    const std::vector<std::pair<int, int>>& getAdjList(int v) const;

    // Get the number of vertices
    int getNumVertices() const;
};

#endif // WEIGHTED_GRAPH_HPP
