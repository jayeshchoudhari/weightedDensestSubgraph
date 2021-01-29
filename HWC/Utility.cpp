#include "Utility.hpp"
#include "Hypergraph.hpp"
#include <cstdio>
#include <iostream> 
#include <vector>
#include <cmath>
using namespace std;

void printVec(const vector<int>& vec){
    for(int x : vec)
        cout << x << " ";
}

void printLnVec(const vector<int>& vec){
    printVec(vec);
    cout << endl;
}

void printVecVec(const vector<vector<int>>& vv){
    for(const vector<int>& v : vv)
        printLnVec(v);
}

void logging(char const* p){
    printf("%s\n", p);
}


template<typename Iter>
pair<unsigned int, unsigned long long> num_nodes_and_edges_induced_tmpl(
        Hypergraph& graph, Iter begin, Iter end) {
    unordered_set<int> nodes(begin, end);
    unordered_set<int> hasEdges;

    for (unordered_set<int>::iterator it = nodes.begin(); it != nodes.end();
            ++it) {

        vector<int> edges;
        graph.edgesContain(*it, &edges);

        for (int edge_id : edges) {
            bool allContained = true;
            vector<int> neighbors;
            graph.edgeVertices(edge_id, &neighbors);

            for (vector<int>::iterator it_neighbor = neighbors.begin();
                    it_neighbor != neighbors.end(); ++it_neighbor) {
                if (nodes.find(*it_neighbor) == nodes.end()) {
                    allContained = false;
                    break;
                }
            }

            if (allContained)
                hasEdges.insert(edge_id);

        }
    }


    return make_pair(nodes.size(), hasEdges.size());
}

pair<unsigned int, unsigned long long> num_nodes_and_edges_induced(
        Hypergraph& graph, vector<int>::iterator begin,
        vector<int>::iterator end) {
    return num_nodes_and_edges_induced_tmpl(graph, begin, end);
}

pair<unsigned int, unsigned long long> num_nodes_and_edges_induced(
        Hypergraph& graph, unordered_set<int>::iterator begin,
        unordered_set<int>::iterator end) {
    return num_nodes_and_edges_induced_tmpl(graph, begin, end);
}

template<typename Iter>
double density_tmpl(Hypergraph& graph, Iter begin, Iter end) {

    pair<int, unsigned long long> induced = num_nodes_and_edges_induced(graph,
            begin, end);

    if (induced.first == 0) {
        return 0;
    }

    return static_cast<double>(induced.second)
            / static_cast<double>(induced.first);
}

double density(Hypergraph& graph, vector<int>::iterator begin,
        vector<int>::iterator end) {
    return density_tmpl(graph, begin, end);
}

double density(Hypergraph& graph, unordered_set<int>::iterator begin,
        unordered_set<int>::iterator end) {
    return density_tmpl(graph, begin, end);
}

int max_iter(const double num_nodes, const double epsilon) {
    return ceil(log(num_nodes) / log1p(epsilon)) + 5; // this is to be sure we have room for rounding error when node count increases
}


// http://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
size_t hash_combine( size_t lhs, size_t rhs ) {
  lhs^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}