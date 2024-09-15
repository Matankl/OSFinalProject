#ifndef MST_STRATEGY_HPP
#define MST_STRATEGY_HPP

#include "WeightedGraph.hpp"
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <memory>
#include <string>
#include <stdexcept>

class MSTStrategy {
public:
    virtual int computeMST(const WeightedGraph& g) = 0;
    virtual ~MSTStrategy() {}
};

class PrimMST : public MSTStrategy {
public:
    int computeMST(const WeightedGraph& g) override;
};

class KruskalMST : public MSTStrategy {
public:
    int computeMST(const WeightedGraph& g) override;

private:
    struct Edge {
        int u, v, weight;
        bool operator<(const Edge& other) const {
            return weight < other.weight;
        }
    };

    std::vector<int> parent, rank;
    int find(int u);
    void unionSets(int u, int v);
};

class MSTFactory {
public:
    static std::unique_ptr<MSTStrategy> createMSTSolver(const std::string& algorithm);
};

#endif // MSTFACTORY_HPP
