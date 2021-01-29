#include "Hypergraph.hpp" 
#include "ApproximateDSFullyDyn.hpp"
#include "GraphScheduler.hpp"
#include "Utility.hpp"
#include "Stats.hpp"
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <string>

// #define DEBUG

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


    GraphScheduler gs(file_name);
    Hypergraph h(true, -1, -1);
    ApproxDSFullyDyn ads(epsilon);
    ads.bindGraph(&h);
    Stats stats(window_sz);
    //stats.bindGraph(&ads);
    std::cout<< "Limit"<< std::numeric_limits<int>::max() <<std::endl;

    while (gs.has_next()) {
        EdgeUpdate edge_up = gs.next_update();
        int edge_id = 0;
        if (!edge_up.is_report) {
            if (edge_up.is_add) {
                edge_id = ads.addEdge(edge_up.vertices);
            } else {
                edge_id = ads.removeEdge(edge_up.vertices);
            }
        }
        assert(edge_id < std::numeric_limits<int>::max() );

        double upper_bound = h.maxEdgeCardinality() * ads.beta() * (1.0 + ads.epsilon());

        if (edge_id != -1){
            stats.exec_op(edge_up.is_add, edge_up.is_report, ads.subgraphSize(),
                    ads.density(), upper_bound, edge_up.timestamp, edge_up.report_label);
        }
    }

    stats.end_op();
    
    return 0;
}
