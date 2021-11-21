#include <iostream>
#include "Graph.h"
#include "WeightedGraph.h"

int main() {
    WeightedGraph weightedGraph({{0, 1, 4}, {0, 7, 8}, {1, 2, 8}, {1, 7, 11}, {2, 3, 7}, {2, 8, 2}, {2, 5, 4}, {3, 4, 9}, {3, 5, 14}, {4, 5, 10}, {5, 6, 2}, {6, 7, 1}, {6, 8, 6}, {7, 8, 7}});
    std::vector<int> mst = weightedGraph.JarnikPrim();
    std::cout << "[Jarnik-Prim] Minimum spanning tree: ";
    for (int i = 0; i < mst.size(); ++i) {
        std::cout << "(" << mst[i] << ", " << i << ") ";
    }
    std::cout << std::endl;

    std::vector<std::pair<int, int>> mstPair = weightedGraph.Kruskal();
    std::cout << "[Kruskal] Minimum spanning tree: ";
    for (int i = 0; i < mstPair.size(); ++i) {
        std::cout << "(" << mstPair[i].first << ", " << mstPair[i].second << ") ";
    }
    std::cout << std::endl;

    weightedGraph = WeightedGraph({{0, 1, 1}, {1, 2, 1}, {2, 3, 2}, {2, 4, 1}, {0, 5, 2}, {5, 6, 3}, {5, 7, 4}});
    std::vector<int> distances = weightedGraph.BellmanFord(0);
    std::cout << "[Bellman-Ford] Distances to other vertices: ";
    for (int i = 0; i < distances.size(); ++i) {
        std::cout << distances[i] << " ";
    }
    std::cout << std::endl;

    distances = weightedGraph.Dijkstra(0);
    std::cout << "[Dijkstra] Distances to other vertices: ";
    for (int i = 0; i < distances.size(); ++i) {
        std::cout << distances[i] << " ";
    }
    std::cout << std::endl;

    weightedGraph = WeightedGraph({{0, 2, -2}, {1, 0, 4}, {1, 2, 3}, {2, 3, 2}, {3, 1, -1}});
    std::vector<std::vector<int>> distanceMatrix = weightedGraph.Johnson();
    std::cout << "[Johnson] Distances from each vertex to each other: " << std::endl;
    for (int i = 0; i < distanceMatrix.size(); ++i) {
        for (int j = 0; j < distanceMatrix.size(); ++j) {
            std::cout << distanceMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    weightedGraph = WeightedGraph({{0, 1, 16}, {0, 2, 13}, {1, 2, 10}, {2, 1, 4}, {1, 3, 12}, {2, 4, 14}, {3, 2, 9}, {4, 3, 7}, {3, 5, 20}, {4, 5, 4}});
    std::cout << "[Ford-Fulkerson] Maximal Flow of Graph: " << weightedGraph.FordFulkerson(0, 5) << std::endl;
    return 0;
}
