#ifndef __UTILITY__
#define __UTILITY__
#include "Hypergraph.hpp"
#include <vector>
#include <unordered_set>
#include <utility>

using namespace std;

void printVec(const  vector<int> &);
void printLnVec(const vector<int> &);
void printVecVec(const  vector<vector<int>> &);
void logging(char const *);


static const int INF_LEVEL = 1000000;

typedef struct DSResult {
    double density;
    std::vector<int> subgraph;
    unsigned long edges_in_subgraph;

    DSResult() {
        density = -1;
        edges_in_subgraph = 0;
    }

    void clear() {
        density = -1;
        edges_in_subgraph = 0;
        subgraph.clear();
    }
} DSResult;

// http://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
size_t hash_combine( size_t lhs, size_t rhs );

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    // return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    return hash_combine(std::hash<T>()(x.first), std::hash<U>()(x.second));
     
  }
};


double density(Hypergraph& graph, vector<int>::iterator subset_begin,
        vector<int>::iterator subset_end);
double density(Hypergraph& graph,
        unordered_set<int>::iterator subset_begin,
        unordered_set<int>::iterator subset_end);

pair<unsigned int, unsigned long long> num_nodes_and_edges_induced(
        Hypergraph& graph, vector<int>::iterator begin,
        vector<int>::iterator end);

pair<unsigned int, unsigned long long> num_nodes_and_edges_induced(
        Hypergraph& graph, unordered_set<int>::iterator begin,
        unordered_set<int>::iterator end);

int max_iter(const double num_nodes, const double epsilon);





#endif