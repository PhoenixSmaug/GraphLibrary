#ifndef GRAPHLIBRARY_WEIGHTEDGRAPH_H
#define GRAPHLIBRARY_WEIGHTEDGRAPH_H

#include <unordered_map>
#include <tuple>
#include "Graph.h"

struct DisjointSet {  // data structure for Kruskal()
private:
    int n;
    std::vector<int> parents, ranks;
public:
    DisjointSet(int n) : n(n) {
        parents = std::vector<int>(n + 1);
        ranks = std::vector<int>(n + 1);
        for (int i = 0; i <= n; ++i) {
            ranks[i] = 0;
            parents[i] = i;
        }
    }

    int find(int v) {
        if (v != parents[v]) {
            parents[v] = find(parents[v]);
        }
        return parents[v];
    }

    void merge(int x, int y) {
        x = find(x), y = find(y);
        if (ranks[x] > ranks[y]) {
            parents[y] = x;
        }
        else {
            parents[x] = y;
        }
        if (ranks[x] == ranks[y]) { ++ranks[y]; }
    }
};

class WeightedGraph : public Graph {
    friend std::ostream &operator<<(std::ostream& os, const WeightedGraph& weightedGraph);
public:
    WeightedGraph(int n);
    WeightedGraph(std::vector<std::pair<int, int>> edges);
    WeightedGraph(std::vector<std::tuple<int, int, float>> edges);

    void addEdge(const int b, const int e);
    void addEdge(const int b, const int e, const float value);
    void addEdge(const std::pair<int, int> edge);
    bool removeEdge(const int b, const int e);
    bool removeEdge(const std::pair<int, int> edge);

    void turnUndirected();

    bool bfs(const int source, const int goal, std::vector<int>& predecessor);

    std::vector<int> JarnikPrim();
    std::vector<std::pair<int, int>> Kruskal();

    std::vector<float> BellmanFord(const int source);
    std::vector<float> Dijkstra(const int source);
    std::vector<int> Dijkstra(const int source, const int goal);
    float DijkstraLength(const int source, const int goal);
    std::vector<std::vector<float>> Johnson();

    float FordFulkerson(const int s, const int t);
    std::vector<int> FordFulkersonPath(const int s, const int t);
private:
    std::unordered_map<std::pair<int,int>, float, pair_hash> weights;
};
std::ostream &operator<<(std::ostream& os, const WeightedGraph& weightedGraph);

#endif //GRAPHLIBRARY_WEIGHTEDGRAPH_H
