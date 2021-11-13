#ifndef GRAPHLIBRARY_GRAPH_H
#define GRAPHLIBRARY_GRAPH_H

#include <vector>
#include <unordered_set>
#include <list>
#include <iostream>

struct pair_hash {  // allow std::pair as key for std::unordered_set
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first * 31+v.second;
    }
};

class Graph {
    friend std::ostream &operator<<(std::ostream& os, const Graph& graph);
public:
    Graph(int n);
    Graph(std::vector<std::pair<int, int>> edges);
    Graph() = default;

    bool existVertex(const int n);  // no input verification for performance reasons
    size_t addVertex();
    bool existEdge(const int b, const int e);
    bool existEdge(const std::pair<int, int> edge);
    void addEdge(const int b, const int e);
    void addEdge(const std::pair<int, int> edge);
    bool removeEdge(const int b, const int e);
    bool removeEdge(const std::pair<int, int> edge);

    void turnUndirected();

    std::vector<bool> bfs(const int source);
    bool bfs(const int source, const int goal);
    std::vector<bool> dfs(const int source);
    bool dfs(const int source, const int goal);
    std::vector<int> topoSort();
protected:
    std::vector<std::list<int>> adjList;
    std::unordered_set<std::pair<int,int>, pair_hash> edges;
};
std::ostream &operator<<(std::ostream& os, const Graph& graph);

#endif //GRAPHLIBRARY_GRAPH_H
