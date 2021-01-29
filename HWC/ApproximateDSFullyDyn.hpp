#ifndef __ApproximateDSFullyDyn__
#define __ApproximateDSFullyDyn__

#include "Hypergraph.hpp"
#include "Utility.hpp"
#include "Decompose.hpp"
#include <utility>
#include <unordered_set>
#include <set>
#include <vector>
#include <assert.h>
#include <stack>
#include <iostream>


class ApproxDSFullyDyn{
public:
    ApproxDSFullyDyn(double eps) : epsilon_(eps) {
        assert(epsilon_ > 0);
        mDegree_ = 0;
        beta_zero_ = 0.1;
        graph_ = NULL; 
    }



    double epsilon(){return epsilon_;}


    void bindGraph(Hypergraph* _hp) {
        graph_ = _hp;
    }

    double beta();
    double density();
    int subgraphSize();

    int addEdge(std::vector<int>& edge);
    int removeEdge(vector<int>& edge);

private:    

    void updateDensity(pair<int, int> pii, double density);
    void checkBetaAlphaSet();

    double beta_zero_ = 0.1;
    double epsilon_;
    std::unordered_map<std::pair<int, int>, Decompose*, pairhash> ab_map_; 
    std::unordered_map<std::pair<int, int>, double, pairhash> density_map_;
    std::set<std::pair<double, std::pair<int, int>>> density_set_;
    int mDegree_;

    Hypergraph * graph_;
};


#endif
