cmake_minimum_required(VERSION 3.20)
project(GraphLibrary)

set(CMAKE_CXX_STANDARD 11)

add_executable(GraphLibrary main.cpp Graph.h Graph.cpp WeightedGraph.cpp WeightedGraph.h)
