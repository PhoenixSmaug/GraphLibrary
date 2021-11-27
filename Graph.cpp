#include <iostream>
#include <stdexcept>
#include <queue>
#include <stack>
#include "Graph.h"

std::ostream &operator<<(std::ostream& os, const Graph& graph) {
    for (int i = 0; i < graph.adjList.size(); ++i) {
        for (auto const& j : graph.adjList[i]) { os << i << " -> " << j << std::endl; }
    }
    return os;
}

Graph::Graph(int n) : adjList(std::vector<std::list<int>>(n)) { }

Graph::Graph(std::vector<std::pair<int, int>> edges) {
    int max = 0;
    for (int i = 0; i < edges.size(); ++i) {
        if (edges[i].first < 0 || edges[i].second < 0) { throw std::invalid_argument("negative value for vertex"); }
        max = (edges[i].first > max) ? edges[i].first : max;
        max = (edges[i].second > max) ? edges[i].second : max;
    }

    this->adjList = std::vector<std::list<int>>(max + 1);
    for (int i = 0; i < edges.size(); ++i) {
        addEdge(edges[i]);
    }
}

bool Graph::existVertex(const int n) { return (n < this->adjList.size()); }

size_t Graph::addVertex() {
    this->adjList.emplace_back(std::list<int>());
    return this->adjList.size() - 1;
}

bool Graph::existEdge(const int b, const int e) {
    return (this->edges.find({b, e}) != this->edges.end());
}

bool Graph::existEdge(const std::pair<int, int> edge) {
    return (this->edges.find(edge) != this->edges.end());
}

void Graph::addEdge(const int b, const int e) {
    this->adjList[b].push_back(e);
    this->edges.insert({b, e});
}

void Graph::addEdge(const std::pair<int, int> edge) {
    this->adjList[edge.first].push_back(edge.second);
    this->edges.insert(edge);
}

bool Graph::removeEdge(const int b, const int e) {
    this->adjList[b].remove(e);
    return this->edges.erase({b, e});
}

bool Graph::removeEdge(const std::pair<int, int> edge) {
    this->adjList[edge.first].remove(edge.second);
    return this->edges.erase(edge);
}

void Graph::turnUndirected() {
    /*
     * Turn Undirected
     * - Complete directed graph to undirected graph
     */

    for (int i = 0; i < this->adjList.size(); ++i) {
        for (auto const& j : this->adjList[i]) {
            if (!this->existEdge(j, i)) {
                addEdge(j, i);
            }
        }
    }
}

std::vector<bool> Graph::bfs(const int source) {
    if (source >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::queue<int> search;
    std::vector<bool> explored(this->adjList.size(), false);
    search.push(source);
    explored[source] = true;

    while (!search.empty()) {
        const int &vertex = search.front();
        search.pop();
        for (auto const& i : this->adjList[vertex]) {
            if (!explored[i]) {
                explored[i] = true;
                search.push(i);
            }
        }
    }

    return explored;
}

bool Graph::bfs(const int source, const int goal) {
    if (source >= this->adjList.size() || goal >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::queue<int> search;
    std::vector<bool> explored(this->adjList.size(), false);
    search.push(source);
    explored[source] = true;

    while (!search.empty()) {
        const int &vertex = search.front();
        search.pop();
        if (vertex == goal) { return true; }
        for (auto const& i : this->adjList[vertex]) {
            if (!explored[i]) {
                explored[i] = true;
                search.push(i);
            }
        }
    }

    return false;
}

std::vector<bool> Graph::dfs(const int source) {
    if (source >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::stack<int> search;
    std::vector<bool> explored(this->adjList.size(), false);
    search.push(source);

    while (!search.empty()) {
        const int &vertex = search.top();
        search.pop();
        if (!explored[vertex]) {
            explored[vertex] = true;
            for (auto const& i : this->adjList[vertex]) {
                search.push(i);
            }
        }
    }

    return explored;
}

bool Graph::dfs(const int source, const int goal) {
    if (source >= this->adjList.size() || goal >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::stack<int> search;
    std::vector<bool> explored(this->adjList.size(), false);
    search.push(source);

    while (!search.empty()) {
        const int &vertex = search.top();
        search.pop();
        if (vertex == goal) { return true; }
        if (!explored[vertex]) {
            explored[vertex] = true;
            for (auto const& i : this->adjList[vertex]) {
                search.push(i);
            }
        }
    }

    return false;
}

std::vector<int> Graph::topoSort() {
    /*
     * Topological Sort
     * - Find topological sort of graph
     * - Returns empty vector if graph contains circles
     * - inspired by https://www.geeksforgeeks.org/topological-sorting-indegree-based-solution
     */

    std::vector<int> indegree(this->adjList.size());
    for (int i = 0; i < this->adjList.size(); ++i) {
        for (auto const& j : this->adjList[i]) {
            ++indegree[j];
        }
    }
    std::queue<int> queue;
    for (int i = 0; i < this->adjList.size(); ++i) {
        if (indegree[i] == 0) {
            queue.push(i);
        }
    }

    int count = 0;
    std::vector<int> order;
    while (!queue.empty()) {
        const int &vertex = queue.front();
        queue.pop();
        order.push_back(vertex);

        for (auto const& i : this->adjList[vertex]) {
            --indegree[i];
            if (indegree[i] == 0) { queue.push(i); }
        }
        ++count;
    }
    if (count != this->adjList.size()) { return std::vector<int>(); }  // graph is not acyclic
    return order;
}