#include "ApproximateDSFullyDyn.hpp"
#include "Utility.hpp"
#include "Hypergraph.hpp"
#include "Decompose.hpp"
#include <vector>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <iomanip> 
#include <stack>
#include <utility>

using namespace std;


int ApproxDSFullyDyn::addEdge(std::vector<int>& edge){

    checkBetaAlphaSet();

    int edge_id = graph_->addEdge(edge);
    if (edge_id == -1)
        return edge_id;


    for (auto it = ab_map_.begin(); it != ab_map_.end(); it++) {
        Decompose* p = it->second;
        double newDensity = p->addEdge(edge, edge_id);
        updateDensity(it->first, newDensity);
    }

    return edge_id; 
}

int ApproxDSFullyDyn::removeEdge(vector<int>& edge){

    checkBetaAlphaSet();


    int edge_id = graph_->removeEdge(edge);
    if (edge_id == -1){
    	// assert(0);
        return edge_id;
    }


    
    for (auto it = ab_map_.begin(); it != ab_map_.end(); it++) {
        Decompose* p = it->second;
        double newDensity = p->removeEdge(edge, edge_id);
        updateDensity(it->first, newDensity);
    }     

    return edge_id;
}

void ApproxDSFullyDyn::checkBetaAlphaSet(){
    if (graph_->maxDegree() + 1> mDegree_) {
        mDegree_ = graph_->maxDegree() + 1;
        while ( beta_zero_ < mDegree_) {
            //TODO FIX
//"Original Line":            const int r = 2;  //graph->maxEdgeCardinality()
            const int r = graph_->maxEdgeCardinality(); // "I changed it- Suman"
            double a = r * (1 + 3 * epsilon_);
            int ab = ceil(a * beta_zero_);
            pair<int, int> pii = make_pair((int)beta_zero_, ab);
            assert (graph_ != NULL);
            if (ab_map_.find(pii) == ab_map_.end() && ab > 1 && (int)beta_zero_ >= 1){
            
                Decompose* p = new Decompose(epsilon_, (int)beta_zero_, ab);
                p->bindGraph(graph_);
                p->newConstruct();
                ab_map_[pii] = p;
            }

            beta_zero_ *= 1 + epsilon_;
        }
    }
    return;
}

void ApproxDSFullyDyn::updateDensity(pair<int, int> pii, double density) {
    if (density_map_.find(pii) == density_map_.end()) {
        density_map_[pii] = density;
        density_set_.insert(make_pair(density, pii));
    } else {
        double oldDensity = density_map_[pii];
        if (density != oldDensity){
            int tmp = density_set_.erase(make_pair(oldDensity, pii));
//            assert(tmp == 1);
            density_set_.insert(make_pair(density, pii));
            density_map_[pii] = density;
        }
    }
    return;
}



double ApproxDSFullyDyn::density() {
    if (ab_map_.empty()) {
        return 0.0;
    }
    return (*density_set_.rbegin()).first;    
}

int ApproxDSFullyDyn::subgraphSize() { 
    if (ab_map_.empty()) {
        return 0;
    }
    return (ab_map_[(*density_set_.rbegin()).second])->subgraphSize();
}

//TODO
double ApproxDSFullyDyn::beta() {
    return 1.1;
}

