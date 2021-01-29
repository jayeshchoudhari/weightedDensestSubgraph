#include "Hypergraph.hpp" 
#include "ApproximateDS.hpp"
#include "Utility.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <ctime>

using namespace std;


int n, m;
// int totalEdgeCardinality;
// double totalEdgeWeight, totalVertexWeight;
// vector<double> wn, we, dn; // weight of vertex, weight of edge, delta of vertex
vector<vector<int>> hyperedge;
// unordered_set<int> Sb;


int main(int argc, char** argv) {

    if (argc != 3) {
        cerr
                << "ERROR Requires 2 parameters. ExecAdd epsilon(double); graph.txt;"
                << endl;
        exit(1);
    }

    double epsilon = atof(argv[1]);
    (void)(freopen(argv[2], "r", stdin) + 1); // remove warnings 

    Hypergraph h(false);
    
    cin >> n >> m;
    assert (n > 0 && m >= 0);

    // totalEdgeCardinality = 0;
    // totalEdgeWeight = 0.0;

    for (int i = 0; i < n; i++) {
        double w;
        cin >> w;
        assert(w > 0);
        // totalVertexWeight += w;
        // wn.push_back(w);
        // dn.push_back(0.0);
    }
    for (int i = 0; i < m; i++) {
        double w;
        cin >> w;
        // we.push_back(w);
        // totalEdgeWeight += w;
        vector<int> edge;
        int ce; //edge cardinality
        cin >> ce;
        assert(ce > 0);
        // totalEdgeCardinality += ce;
        for (int j = 0; j < ce; j++) {
            int vid;
            cin >> vid;
            assert(vid >= 0 && vid < n);
            edge.push_back(vid);
            // dn[vid] += (1.0 * w) / ce;
        }
        //hyperedge.push_back(edge);
        h.addEdge(edge);
    }


    
    std::clock_t start = std::clock();

    ApproxDS ads(epsilon);
    ads.bindGraph(&h);

    // for (auto e : hyperedge)
        // ads.addEdge(e);

    ads.Construct();

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Duration: " << 1000 * duration << " ";

    // double upper_bound = h.maxEdgeCardinality() * ads.beta() * (1.0 + ads.epsilon());

    cout << "Density: " << ads.subgraphDensity() << " VertexSize: " << ads.subgraphSize()  << " EdgeSize: " << ads.subgraphEdgeSize() << "\n";

    fclose (stdin);

    return 0;
}