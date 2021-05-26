#ifndef __ApproximateDS__
#define __ApproximateDS__

#include "Hypergraph.hpp"
#include "Utility.hpp"
#include <unordered_set>
#include <vector>
#include <assert.h>



class ApproxDS{
public:

    ApproxDS(double eps) { 
        epsilon_ = eps;
        edges_densest_subgraph_ = 0;
        beta_ = 0.25/ (1.0 + eps); 
        assert(epsilon_ > 0); 
    }



    double subgraphDensity() const {
        if (densest_subgraph_set_.size() == 0) {
            return 0.0;
        }
        return static_cast<double>(edges_densest_subgraph_)
                / static_cast<double>(densest_subgraph_set_.size());
    }    

    unsigned int subgraphSize() const {
        return densest_subgraph_set_.size();
    }

    int subgraphEdgeSize() const {
        return edges_densest_subgraph_;
    }

    std::vector<int> subgraph(){
        std::vector<int> v;
        for (int x : densest_subgraph_set_)
            v.push_back(x);
        return v;
    }



    void bindGraph(Hypergraph* _hp) {
        graph_ = _hp;
    }

    int addEdge(std::vector<int>& edge);

    void Construct();

    void removeEdge(const int edge_id){ 
        //graph_->removeEdge(edge_id);
        if (edge_id > 0){
            ; // do noting, eliminating compiling warnings
        }
        assert(0);
    }

    double beta(){return beta_;}
    double epsilon(){return epsilon_;}

private:



    double beta_, epsilon_;

    std::unordered_set<int> densest_subgraph_set_;
    int edges_densest_subgraph_;

    std::unordered_map<int, int> level_map_;
    std::unordered_map<int, int> level_degree_map_;

    Hypergraph * graph_;
};


#endif