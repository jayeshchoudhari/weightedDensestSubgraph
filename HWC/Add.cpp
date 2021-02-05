#include "Hypergraph.hpp" 
#include "ApproximateDS.hpp"
#include "GraphScheduler.hpp"
#include "Utility.hpp"
#include "Stats.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>


using namespace std;



int main(int argc, char** argv) {

    if (argc != 4) {
        cerr
                << "ERROR Requires 3 parameters. ExecAdd epsilon(double); graph-udates.txt; stat_window_size;"
                << endl;
        exit(1);
    }

    double epsilon = atof(argv[1]);
    string file_name(argv[2]);
    const int window_sz = atoi(argv[3]);
    // the r parameter is modified to have the rank of the
    // hypergraph as an input (It can be 2 or the actual rank)
    const  int r = atoi(argv[3]);


    GraphScheduler gs(file_name);
    Hypergraph h(false, r);
    ApproxDS ads(epsilon);
    ads.bindGraph(&h);
    Stats stats(window_sz);
    
    while (gs.has_next()) {
        EdgeUpdate edge_up = gs.next_update();
        int edge_id = 0;
        if (!edge_up.is_report) {
            if (edge_up.is_add) {
                edge_id = ads.addEdge(edge_up.vertices);
            } else {
                continue;
            }
        }
        assert(edge_id < std::numeric_limits<int>::max() );

        double upper_bound = h.maxEdgeCardinality() * ads.beta() * (1.0 + ads.epsilon());
        if (edge_id != -1){
            stats.exec_op(edge_up.is_add, edge_up.is_report, ads.subgraphSize(),
                          ads.subgraphDensity(), upper_bound, edge_up.timestamp, edge_up.report_label);
        }

//        if (edge_id != -1)
//            stats.exec_op(edge_up.is_add, ads.subgraphSize(), ads.subgraphDensity(), upper_bound, edge_up.timestamp);
    }

    stats.end_op();
    
    return 0;
}