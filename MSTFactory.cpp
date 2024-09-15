#include "MSTFactory.hpp"

// Implementation of Prim's MST algorithm
// This function calculates the Minimum Spanning Tree (MST) for the graph 'g' using Prim's algorithm
// It returns the total weight of the MST
int PrimMST::computeMST(const WeightedGraph& g) {
    int n = g.getNumVertices();  // Get the number of vertices in the graph

    // Key array to store the minimum weight to reach each vertex, initialized to infinity
    std::vector<int> key(n + 1, std::numeric_limits<int>::max());

    // Boolean array to keep track of vertices included in the MST
    std::vector<bool> inMST(n + 1, false);

    // Set the key for the first vertex to 0, so it's picked first
    key[1] = 0;

    // Min-priority queue to pick the vertex with the smallest key value (pair of <key, vertex>)
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    pq.push({0, 1});  // Push the first vertex with key 0

    int totalWeight = 0;  // Initialize the total weight of the MST

    // While there are vertices not yet included in the MST
    while (!pq.empty()) {
        int u = pq.top().second;  // Get the vertex with the smallest key value
        pq.pop();

        // If this vertex is already in the MST, skip it
        if (inMST[u]) continue;

        // Include vertex 'u' in the MST and add its key value to the total weight
        inMST[u] = true;
        totalWeight += key[u];

        // Iterate over all adjacent vertices of 'u'
        for (const auto &[v, weight] : g.getAdjList(u)) {
            // If 'v' is not in the MST and the weight of the edge (u, v) is smaller than the current key for 'v'
            if (!inMST[v] && weight < key[v]) {
                key[v] = weight;  // Update the key with the smaller weight
                pq.push({key[v], v});  // Push the updated key and vertex into the priority queue
            }
        }
    }
    return totalWeight;  // Return the total weight of the MST
}

// Find method for Kruskal's MST
// This function finds the representative of the set containing vertex 'u' using path compression
int KruskalMST::find(int u) {
    // If 'u' is not its own parent, recursively find the root parent and compress the path
    if (parent[u] != u) parent[u] = find(parent[u]);
    return parent[u];
}

// UnionSets method for Kruskal's MST
// This function unites two disjoint sets containing 'u' and 'v' by rank
void KruskalMST::unionSets(int u, int v) {
    u = find(u);  // Find the root of set containing 'u'
    v = find(v);  // Find the root of set containing 'v'

    // Union by rank: attach the tree with lower rank under the tree with higher rank
    if (rank[u] < rank[v]) parent[u] = v;
    else if (rank[u] > rank[v]) parent[v] = u;
    else parent[v] = u, rank[u]++;  // If ranks are equal, increase the rank of 'u'
}

// Implementation of Kruskal's MST algorithm
// This function calculates the MST for the graph 'g' using Kruskal's algorithm
// It returns the total weight of the MST
int KruskalMST::computeMST(const WeightedGraph& g) {
    int n = g.getNumVertices();  // Get the number of vertices in the graph

    std::vector<Edge> edges;  // Vector to store all the edges of the graph
    parent.resize(n + 1);     // Resize the parent array for union-find operations
    rank.assign(n + 1, 0);    // Initialize the rank array with zeros

    // Iterate over all vertices and their adjacency lists to gather all edges
    for (int u = 1; u <= n; ++u) {
        for (const auto &[v, w] : g.getAdjList(u)) {
            if (u < v) edges.push_back({u, v, w});  // Add edge (u, v) with weight 'w'
        }
    }

    // Sort all edges by weight in non-decreasing order
    std::sort(edges.begin(), edges.end());

    // Initialize each vertex to be its own parent in the union-find structure
    for (int i = 1; i <= n; ++i) parent[i] = i;

    int totalWeight = 0;  // Initialize the total weight of the MST

    // Iterate over all sorted edges
    for (const auto &edge : edges) {
        // If the two vertices of the edge belong to different sets, add the edge to the MST
        if (find(edge.u) != find(edge.v)) {
            totalWeight += edge.weight;  // Add the edge's weight to the total weight
            unionSets(edge.u, edge.v);   // Union the sets of the two vertices
        }
    }

    return totalWeight;  // Return the total weight of the MST
}

// Factory method to create MST solver based on the algorithm name
// This function returns an instance of either PrimMST or KruskalMST based on the provided algorithm name
// Throws an exception if an unsupported algorithm is requested
std::unique_ptr<MSTStrategy> MSTFactory::createMSTSolver(const std::string& algorithm) {
    if (algorithm == "Prim") {
        // If the algorithm is "Prim", return an instance of PrimMST
        return std::make_unique<PrimMST>();
    } else if (algorithm == "Kruskal") {
        // If the algorithm is "Kruskal", return an instance of KruskalMST
        return std::make_unique<KruskalMST>();
    } else {
        // If the algorithm is not recognized, throw an exception
        throw std::invalid_argument("Unsupported MST algorithm");
    }
}
