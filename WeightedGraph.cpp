#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include "WeightedGraph.h"

std::ostream &operator<<(std::ostream& os, const WeightedGraph& weightedGraph) {
    for (int i = 0; i < weightedGraph.adjList.size(); ++i) {
        for (auto const& j : weightedGraph.adjList[i]) { std::cout << i << " -> " << j << " : " << weightedGraph.weights.at({i, j}) << std::endl; }
    }
    return os;
}

WeightedGraph::WeightedGraph(int n) : Graph(n) {}

WeightedGraph::WeightedGraph(std::vector<std::pair<int, int>> edges) {
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

WeightedGraph::WeightedGraph(std::vector<std::tuple<int, int, int>> edges) {
    int max = 0;
    for (int i = 0; i < edges.size(); ++i) {
        if (std::get<0>(edges[i]) < 0 || std::get<1>(edges[i]) < 0) { throw std::invalid_argument("negative value for vertex"); }
        max = (std::get<0>(edges[i]) > max) ? std::get<0>(edges[i]) : max;
        max = (std::get<1>(edges[i]) > max) ? std::get<1>(edges[i]) : max;
    }

    this->adjList = std::vector<std::list<int>>(max + 1);
    for (int i = 0; i < edges.size(); ++i) {
        addEdge(std::get<0>(edges[i]), std::get<1>(edges[i]), std::get<2>(edges[i]));
    }
}

void WeightedGraph::addEdge(const int b, const int e) {
    Graph::addEdge(b, e);
    weights.insert({{b, e}, 1});
}

void WeightedGraph::addEdge(const std::pair<int, int> edge) {
    Graph::addEdge(edge.first, edge.second);
    weights.insert({edge, 1});
}

void WeightedGraph::addEdge(const int b, const int e, const int value) {
    Graph::addEdge(b, e);
    weights.insert({{b, e}, value});
}

bool WeightedGraph::removeEdge(const int b, const int e) {
    weights.erase({b, e});
    return Graph::removeEdge(b, e);
}

bool WeightedGraph::removeEdge(const std::pair<int, int> edge) {
    weights.erase(edge);
    return Graph::removeEdge(edge);
}

void WeightedGraph::turnUndirected() {
    for (int i = 0; i < this->adjList.size(); ++i) {
        for (auto const& j : this->adjList[i]) {
            if (!this->existEdge(j, i)) {
                addEdge(j, i, this->weights.at({i, j}));
            }
        }
    }
}

bool WeightedGraph::bfs(const int source, const int goal, std::vector<int>& predecessor) {
    if (source >= this->adjList.size() || goal >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::queue<int> search;
    std::vector<bool> explored(this->adjList.size(), false);
    search.push(source);
    explored[source] = true;
    predecessor[source] = -1;

    while (!search.empty()) {
        const int &vertex = search.front();
        search.pop();
        for (auto const& i : this->adjList[vertex]) {
            if (!explored[i] && (this->weights.at({vertex, i}) > 0)) {
                if (i == goal) {
                    predecessor[i] = vertex;
                    return true;
                }

                explored[i] = true;
                search.push(i);
                predecessor[i] = vertex;
            }
        }
    }

    return false;
}

std::vector<int> WeightedGraph::JarnikPrim() {
    /*
     * Jarnik-Prim Algorithm
     * - Find minimal spanning tree
     * - Use with undirected graphs, no input validation
     * - inspired by https://www.geeksforgeeks.org/prims-algorithm-using-priority_queue-stl/
     */

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> priorityQueue;
    std::vector<int> keys(this->adjList.size(), std::numeric_limits<int>::max());
    std::vector<int> parents(this->adjList.size(), -1);
    std::vector<int> inMST(this->adjList.size(), false);

    priorityQueue.push({0, 0});
    keys[0] = 0;

    while (!priorityQueue.empty()) {
        int min = priorityQueue.top().second;  // take greedy edge with minimal weight
        priorityQueue.pop();
        if (inMST[min]) { continue; }
        inMST[min] = true;

        for (auto const& i : this->adjList[min]) {
            int weight = this->weights.at({min, i});
            if (!inMST[i] && keys[i] > weight) {
                keys[i] = weight;
                priorityQueue.push({keys[i], i});
                parents[i] = min;
            }
        }
    }

    parents.erase(parents.begin());
    return parents;  // MST with edges (parents[i], i)
}

std::vector<std::pair<int, int>> WeightedGraph::Kruskal() {
    /*
     * Kruskal Algorithm
     * - Find minimal spanning tree
     * - Use with undirected graphs, no input validation
     * - inspired by https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
     */

    std::vector<std::pair<int, int>> collector;
    std::vector<std::pair<int, std::pair<int, int>>> edges;  // {value, {v, w}} with v < w
    for (int i = 0; i < this->adjList.size(); ++i) {
        for (auto const& j : this->adjList[i]) {
            if (i < j) { edges.push_back({this->weights.at({i, j}), {i, j}}); }
        }
    }
    std::sort(edges.begin(), edges.end());
    DisjointSet disjointSet(this->adjList.size());

    for (auto const& i : edges) {
        int v = i.second.first;
        int w = i.second.second;
        int setV = disjointSet.find(v);
        int setW = disjointSet.find(w);

        if (setV != setW) {  // guarantee acyclic graph
            collector.push_back({v, w});
            disjointSet.merge(setV, setW);
        }
    }

    return collector;
}

std::vector<int> WeightedGraph::BellmanFord(const int source) {
    /*
     * Bellman-Ford Algorithm
     * - Find distance to each vertex in graph
     * - Returns empty vector if graph is not conservative (contains cycles with negative length)
     * - Not reachable vertices are returned as MAX_INT
     */

    if (source >= this->adjList.size()) { throw std::invalid_argument("node not in graph"); }
    std::vector<int> distances(this->adjList.size(), std::numeric_limits<int>::max());
    distances[source] = 0;

    for (int i = 0; i < this->adjList.size() - 1; ++i) {
        for (int j = 0; j < this->adjList.size(); ++j) {
            for (auto const& k : this->adjList[j]) {
                if (distances[j] != std::numeric_limits<int>::max() && distances[j] + this->weights.at({j, k}) < distances[k]) {
                    distances[k] = distances[j] + this->weights.at({j, k});
                }
            }
        }
    }

    for (int i = 0; i < this->adjList.size(); ++i) {  // graph is not conservative
        for (auto const& j : this->adjList[i]) {
            if (distances[i] != std::numeric_limits<int>::max() && distances[i] + this->weights.at({i, j}) < distances[j]) {
                return std::vector<int>();
            }
        }
    }

    return distances;
}

std::vector<int> WeightedGraph::Dijkstra(const int source) {
    /*
     * Dijkstra Algorithm
     * - Find distance to each vertex in graph
     * - Only positive weights allowed, no input validation
     * - Not reachable vertices are returned as MAX_INT
     * - inspired by https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/
     */

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> priorityQueue;
    std::vector<int> distances(this->adjList.size(), std::numeric_limits<int>::max());
    priorityQueue.push({0, source});
    distances[source] = 0;

    while (!priorityQueue.empty()) {
        int min = priorityQueue.top().second;
        priorityQueue.pop();

        for (auto const& i : this->adjList[min]) {
            int weight = this->weights.at({min, i});
            if (distances[i] > distances[min] + weight) {
                distances[i] = distances[min] + weight;
                priorityQueue.push({distances[i], i});
            }
        }
    }

    return distances;
}

std::vector<int> WeightedGraph::Dijkstra(const int source, const int goal) {
    /*
     * Dijkstra Algorithm
     * - Find the shortest path from source to goal
     * - Only positive weights allowed, no input validation
     * - Returns empty vector if goal is not reachable
     */

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> priorityQueue;
    std::vector<int> distances(this->adjList.size(), std::numeric_limits<int>::max());
    std::vector<int> predecessors(this->adjList.size(), -1);
    priorityQueue.push({0, source});
    distances[source] = 0;

    while (!priorityQueue.empty()) {
        int min = priorityQueue.top().second;
        priorityQueue.pop();

        if (min == goal) {
            std::vector<int> path;
            if (predecessors[min] == -1 || min == source) { return std::vector<int>(); }
            while (min != -1) {
                path.push_back(min);
                min = predecessors[min];
            }
            return path;
        }

        for (auto const &i: this->adjList[min]) {
            int weight = this->weights.at({min, i});
            if (distances[i] > distances[min] + weight) {
                distances[i] = distances[min] + weight;
                priorityQueue.push({distances[i], i});
                predecessors[i] = min;
            }
        }
    }

    return std::vector<int>();
}

int WeightedGraph::DijkstraLength(const int source, const int goal) {
    /*
     * Dijkstra Algorithm
     * - Find the shortest distance from source to goal
     * - Only positive weights allowed, no input validation
     * - Returns 0 if goal is not reachable
     */

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> priorityQueue;
    std::vector<int> distances(this->adjList.size(), std::numeric_limits<int>::max());
    std::vector<int> predecessors(this->adjList.size(), -1);
    priorityQueue.push({0, source});
    distances[source] = 0;

    while (!priorityQueue.empty()) {
        int min = priorityQueue.top().second;
        priorityQueue.pop();

        if (min == goal) {
            int path = 0;
            if (predecessors[min] == -1 || min == source) { return 0; }
            while (predecessors[min] != -1) {
                path += this->weights.at({predecessors[min], min});
                min = predecessors[min];
            }
            return path;
        }

        for (auto const &i: this->adjList[min]) {
            int weight = this->weights.at({min, i});
            if (distances[i] > distances[min] + weight) {
                distances[i] = distances[min] + weight;
                priorityQueue.push({distances[i], i});
                predecessors[i] = min;
            }
        }
    }

    return 0;
}

std::vector<std::vector<int>> WeightedGraph::Johnson() {
    /*
     * Johnson's Algorithm
     * - Find shortest distance between all pair of vertices in graph
     * - Returns empty 2d vector if graph is not conservative (contains cycles with negative length)
     * - Not reachable vertices are returned as MAX_INT
     */

    WeightedGraph weightedGraph = *this;
    weightedGraph.addVertex();
    for (int i = 0; i < weightedGraph.adjList.size() - 1; ++i) {
        weightedGraph.adjList[weightedGraph.adjList.size() - 1].push_back(i);
        weightedGraph.weights.insert({{weightedGraph.adjList.size() - 1, i}, 0});
    }

    std::vector<int> h = weightedGraph.BellmanFord(weightedGraph.adjList.size() - 1);
    if (h.empty()) { return std::vector<std::vector<int>>(weightedGraph.adjList.size() - 1, std::vector<int>()); }

    for (int i = 0; i < weightedGraph.adjList.size(); ++i) {
        for (auto const &j: weightedGraph.adjList[i]) {
            weightedGraph.weights.at({i, j}) = weightedGraph.weights.at({i, j}) + h[i] - h[j];
        }
    }

    weightedGraph.adjList.pop_back();
    std::vector<std::vector<int>> distances(weightedGraph.adjList.size());
    for (int i = 0; i < weightedGraph.adjList.size(); ++i) {
        distances[i] = weightedGraph.Dijkstra(i);
    }

    for (int i = 0; i < weightedGraph.adjList.size(); ++i) {
        for (int j = 0; j < weightedGraph.adjList.size(); ++j) {
            if (distances[i][j] != std::numeric_limits<int>::max()) {
                distances[i][j] += h[j] - h[i];
            }
        }
    }

    return distances;
}

int WeightedGraph::FordFulkerson(const int s, const int t) {
    /*
     * Ford-Fulkerson Algorithm
     * - Find maximal flow between source (s) and sink (t)
     * - Return 0 if no flow exists
     * - inspired by https://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
     */

    WeightedGraph residualGraph = *this;
    for (int i = 0; i < this->adjList.size(); ++i) {
        for (int j = 0; j < this->adjList.size(); ++j) {
            if (residualGraph.edges.count({i, j}) == 0) {
                residualGraph.addEdge(i, j, 0);
            }
        }
    }

    std::vector<int> predecessors(this->adjList.size());
    int maxFlow = 0;

    while(residualGraph.bfs(s, t, predecessors)) {
        int pathFlow = std::numeric_limits<int>::max();
        for (int i = t; i != s; i = predecessors[i]) {
            int j = predecessors[i];
            pathFlow = std::min(pathFlow, residualGraph.weights.at({j, i}));
        }

        for (int i = t; i != s; i = predecessors[i]) {
            int j = predecessors[i];
            residualGraph.weights.at({j, i}) -= pathFlow;
            residualGraph.weights.at({i, j}) += pathFlow;
        }
        maxFlow += pathFlow;
    }

    return maxFlow;
}

std::vector<int> WeightedGraph::FordFulkersonPath(const int s, const int t) {
    /*
    * Ford-Fulkerson Algorithm
    * - Find used edges of maximal flow between source (s) and sink (t)
    * - Return empty vector if no flow exists
    * - inspired by https://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
    */

    WeightedGraph residualGraph = *this;
    for (int i = 0; i < this->adjList.size(); ++i) {
        for (int j = 0; j < this->adjList.size(); ++j) {
            if (residualGraph.edges.count({i, j}) == 0) {
                residualGraph.addEdge(i, j, 0);
            }
        }
    }

    std::vector<int> predecessors(this->adjList.size());
    std::set<int> path;
    path.insert(s);

    while(residualGraph.bfs(s, t, predecessors)) {
        int pathFlow = std::numeric_limits<int>::max();
        for (int i = t; i != s; i = predecessors[i]) {
            int j = predecessors[i];
            pathFlow = std::min(pathFlow, residualGraph.weights.at({j, i}));
        }

        for (int i = t; i != s; i = predecessors[i]) {
            int j = predecessors[i];
            residualGraph.weights.at({j, i}) -= pathFlow;
            residualGraph.weights.at({i, j}) += pathFlow;
            path.insert(i);
        }
    }

    return std::vector<int>(path.begin(), path.end());
}