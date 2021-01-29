#ifndef __DECOMPOSE__
#define __DECOMPOSE__

#include "Hypergraph.hpp"
#include "ApproximateDSFullyDyn.hpp"
#include "Utility.hpp"
#include <utility>
#include <unordered_set>
#include <set>
#include <vector>
#include <assert.h>
#include <stack>
#include <iostream>


class Decompose{
public:
    Decompose(double e, int b, int ab) : epsilon_(e), beta_(b), alpha_beta_(ab), graph_(NULL){}

    double addEdge(vector<int>&, int);
    double removeEdge(vector<int>&, int);

    void bindGraph(Hypergraph* _hp) {
        graph_ = _hp;
    }
    
    void newConstruct();

    double density();

    int subgraphSize();

private:
    void Promote(int, stack<int>&);

    void Demote(int, stack<int>&);

    void Check(std::vector<int>&, int);

    void Add(int u, int l, int edge_id);
    
    void Erase(int u, int l, int edge_id);

    void updateTao();


    void updateDensity();

    double epsilon_;

    int beta_, alpha_beta_;

    int tao_;

    unordered_set<int> touchedLevel;

    // l(u)
    std::unordered_map<int, int> level_map_;  
    //E_u(l)
    std::unordered_map<std::pair<int, int>, std::vector<int>, pairhash> rmEdge_map_;

    std::unordered_map<int, int> level_vertices_, level_edges_;

    std::unordered_map<int, double> level_density_;
    std::set<std::pair<double, int>> level_density_set_;

    Hypergraph * graph_;
};




#endif
