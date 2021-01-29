#include "ApproximateDS.hpp"
#include "Utility.hpp"
#include "Hypergraph.hpp"

#include <vector>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include <iostream>
#include <iomanip> 
#include <stack>

using namespace std;

void Find(const Hypergraph* graph, const double beta, const double epsilon,
        DSResult * result, unordered_map<int, int>* level_map, bool* valid_orientation,
        unordered_map<int, int>* level_degree_map) {


    assert(level_degree_map != NULL);
    level_degree_map->clear();
    

    vector<int> nodes;
    graph->vertices(&nodes);
    unordered_set<int> present(nodes.begin(), nodes.end());

    unordered_set<int> to_delete;
    unordered_set<int> to_delete_next;

    Hypergraph m_graph = *graph;

    result->density = m_graph.graphDensity();
    m_graph.vertices(&result->subgraph);
    level_map->clear();

    //cout << " [IN FIND] begin dens " << result->density << endl;

    const int max_iter_num = max_iter(nodes.size(), epsilon);

    const int max_r = graph->maxEdgeCardinality();

    const double threshold = 1.0 * max_r * beta * (1.0 + epsilon);
    int cutDegree = (int)ceil(threshold);
    
    assert(threshold > 0.0);

    vector<int> v_neighbors;

    for (vector<int>::iterator nodes_it = nodes.begin();
            nodes_it != nodes.end(); ++nodes_it) {
        int v = *nodes_it;
        if (m_graph.degree(v) < cutDegree) {
            to_delete.insert(v);
        }
    }



    for (int iter = 1; iter <= max_iter_num; ++iter) {
        if (to_delete.empty()) {
            break;
        }

        to_delete_next.clear();

        for (unordered_set<int>::iterator to_delete_it = to_delete.begin();
                to_delete_it != to_delete.end(); ++to_delete_it) {
            int v = *to_delete_it;
           (*level_degree_map)[v] = m_graph.degree(v);
        }
        
        for (unordered_set<int>::iterator to_delete_it = to_delete.begin();
                to_delete_it != to_delete.end(); ++to_delete_it) {
            int v = *to_delete_it;

            if (present.find(v) == present.end()) {
                continue;
            }

            present.erase(v);
            v_neighbors.clear();
            m_graph.neighborsOf(v, &v_neighbors);
            (*level_map)[v] = iter;
            m_graph.removeVertex(v);

            for (vector<int>::iterator neighbor_v_it = v_neighbors.begin();
                    neighbor_v_it != v_neighbors.end(); ++neighbor_v_it) {
                int neighbor_v = *neighbor_v_it;

                if (present.find(neighbor_v) != present.end() && to_delete.find(neighbor_v) == to_delete.end()) {
                    if (m_graph.degree(neighbor_v) < cutDegree) {
                        to_delete_next.insert(neighbor_v);
                    }
                }
            }
        }

        double new_density = m_graph.graphDensity();
        if (new_density > result->density) {
            result->density = new_density;
            result->subgraph.clear();
            m_graph.vertices(&result->subgraph);
            //cout << " [IN FIND] improved density " << result->density << endl;
            //cout << " [IN FIND] SIZE S" << result->subgraph.size() << endl;
            //cout << " [IN FIND] EDGES in curr subgraph" << m_graph.num_edges()
            //      << endl;
            //cout << " [IN FIND] NODES in curr subgraph" << m_graph.num_nodes()
            //      << endl;

        }

        to_delete.clear();
        to_delete.insert(to_delete_next.begin(), to_delete_next.end());
    }

    *valid_orientation = present.empty();
}

// Starting from beta_zero, runs Find until it fails.
// If provided it will use the current_solution graph as starting point to set beta.
template<typename Iter>
void Incremental_tmpl(Hypergraph * graph, const double beta_zero,
        const double epsilon, const Iter& current_subgraph_begin,
        const Iter& current_subgraph_end, DSResult * result,
        unordered_map<int, int>* level_map,
        unordered_map<int, int>* level_degree_map) {

    DSResult optimal;
    double beta = beta_zero;

    double curr_dens = density(*graph, current_subgraph_begin,
            current_subgraph_end);
    //cout << "fd " << beta_zero <<curr_dens << endl;
    assert(beta_zero > 0 || curr_dens > 0);

    if (curr_dens > 0) {
        optimal.subgraph.clear();
        optimal.subgraph.assign(current_subgraph_begin, current_subgraph_end);
        optimal.density = curr_dens;
    } else {
        optimal.subgraph.clear();
        optimal.density = -1;
    }

    beta = max(beta_zero, curr_dens * (1.0 + epsilon));

    while (true) {
        bool valid_orientation = false;
        DSResult new_res;
        //cerr << "Find: " << beta << endl;
        //cerr << "Curr dens: " << optimal.density << endl;
        //cerr << "SIZE ds: " << optimal.subgraph.size() << endl;
        //cerr << "Edges graphs: " << graph.num_edges() << endl;

        Find(graph, beta, epsilon, &new_res, level_map,
                &valid_orientation, level_degree_map);
        if (new_res.density >= beta && new_res.density > optimal.density) {
            optimal.density = new_res.density;
            //cout << " [IN ITERATIVE] " << new_res.density << endl;
            optimal.subgraph.clear();
            optimal.subgraph.assign(new_res.subgraph.begin(),
                    new_res.subgraph.end());
            beta = new_res.density * (1.0 + epsilon);
        } else {
             assert(valid_orientation); //this must be true for the theory.    
             assert(level_degree_map->size() == graph->numOfVertices());
            
            break;
        }
    }

    assert(optimal.density >= 0);
    assert(optimal.density * (1 + epsilon) >= beta);
    *result = optimal;

    assert(graph->numOfEdges() == 0 || optimal.density > 0);

    //cout << "END INCREMENTAL: " << endl;
    //cout << "Curr dens: " << result->density << endl;
    //cout << "SIZE ds: " << result->subgraph.size() << endl;
}

void Incremental(Hypergraph * graph, const double beta_zero,
        const double epsilon, vector<int>::iterator current_subgraph_begin,
        vector<int>::iterator current_subgraph_end, DSResult * result,
        unordered_map<int, int>* level_map,
        unordered_map<int, int>* level_degree_map) {
    Incremental_tmpl(graph, beta_zero, epsilon, current_subgraph_begin,
            current_subgraph_end, result, level_map,
            level_degree_map);
}


void Incremental(Hypergraph * graph, const double beta_zero,
        const double epsilon, unordered_set<int>::iterator current_subgraph_begin,
        unordered_set<int>::iterator current_subgraph_end, DSResult * result,
        unordered_map<int, int>* level_map,
        unordered_map<int, int>* level_degree_map) {
    Incremental_tmpl(graph, beta_zero, epsilon, current_subgraph_begin,
            current_subgraph_end, result, level_map,
            level_degree_map);
}



void ApproxDS::Construct() {
    DSResult result;
    Incremental(graph_, beta_, epsilon_, densest_subgraph_set_.begin(),
            densest_subgraph_set_.end(), &result, &level_map_, &level_degree_map_);

    densest_subgraph_set_.clear();
    densest_subgraph_set_.insert(result.subgraph.begin(),
            result.subgraph.end());

    pair<int, unsigned long long> nodes_edges_pair =
            num_nodes_and_edges_induced(*graph_, densest_subgraph_set_.begin(),
                    densest_subgraph_set_.end());
    edges_densest_subgraph_ = nodes_edges_pair.second;

    beta_ = result.density * (1 + epsilon_);
}




int ApproxDS::addEdge(vector<int>& edge) {

    
    int edge_id = graph_->addEdge(edge);

    if (edge_id == -1)
        return -1;


    bool allContained = true;
    for (const int u : edge) {
        if (densest_subgraph_set_.find(u) == densest_subgraph_set_.end()) {
            allContained = false;
            break;
        }
    }
    if (allContained)
        ++edges_densest_subgraph_;


    int min_level = INF_LEVEL;
    
    bool hasNewNodes = false;
    for (const int u : edge) {
        if (level_map_.find(u) == level_map_.end()) {
            hasNewNodes = true;
            break;
        }
    }

    
    if (!hasNewNodes) {
        for (const int u : edge) 
            min_level = min(min_level, level_map_[u]);
    } else {
        min_level = 1;
    }

    assert(min_level != INF_LEVEL);

    for (const int u : edge){
        // For new comers, degree at most 1
        if (level_map_.find(u) == level_map_.end()) {
            level_map_[u] = 1;
            level_degree_map_[u] = 0;
        } 

        if(level_map_[u] == min_level){
            level_degree_map_[u] += 1;
        }
    }


    const int max_r = graph_->maxEdgeCardinality();

    const double threshold = 1.0 * max_r * beta_ * (1.0 + epsilon_);

    int cutDegree = ceil(threshold);
    int max_iter_num = max_iter(graph_->numOfVertices(), epsilon_);

    stack<int> badVertices;
    for (int x : edge)
        badVertices.push(x);

    //cerr << " hi " << endl;

    while(!badVertices.empty()) {
        int u = badVertices.top();
        badVertices.pop();
        int curLevel = level_map_[u];

        while(level_degree_map_[u] >= cutDegree && curLevel < max_iter_num) {    
            curLevel += 1;
            vector<int> edges;
            graph_->edgesContain(u, &edges);
            int cnt = 0;
            for (int edge_id : edges) {
                vector<int> vertices;
                bool allAbove = true;
                graph_->edgeVertices(edge_id, &vertices);
                for (int v : vertices) {
                    if (u == v)
                        continue;
                    if (level_map_[v] < curLevel){
                        allAbove = false;
                        break;
                    }
                }
                if(allAbove)
                    ++cnt;
            }
            level_degree_map_[u] = cnt;
        }

        if (level_map_[u] >= max_iter_num){
            // cerr << " recontruct " << beta_ << " threshold " << threshold << " cut "  << cutDegree << endl;
            // cerr << " node u " << u << " level " << level_map_[u] << " degree "<<level_degree_map_[u] << endl;          
            Construct();
            return edge_id;
        }

        if (curLevel > level_map_[u]) {
            vector<int> edges;
            graph_->edgesContain(u, &edges);
            for (int edge_id : edges) {
                int min_level = INF_LEVEL;
                vector<int> vertices;
                graph_->edgeVertices(edge_id, &vertices);

                // singleton hyperedge, does not change anything
                if (vertices.size() == 1)
                    continue;

                for (int v : vertices) {
                    if (u != v)
                        min_level = min(min_level, level_map_[v]);
                } 
                assert(min_level != INF_LEVEL && min_level > 0);
                if (min_level <= level_map_[u]){
                    // other nodes affect, do nothing
                    ;
                } else {
                    for (int v : vertices){
                        if (u != v) {
                            assert (level_map_[v] >= min_level);
                            if (level_map_[v] == min_level){
                                level_degree_map_[v]++;
                                badVertices.push(v);
                            }
                        }
                    }
                }

            }

            level_map_[u] = curLevel;

            {
                // recalculate b[u]
                vector<int> edges;
                graph_->edgesContain(u, &edges);
                int cnt = 0;
                for (int edge_id : edges) {
                    vector<int> vertices;
                    graph_->edgeVertices(edge_id, &vertices);
                    bool allAbove = true;
                    for (int v : vertices) {
                        if (level_map_[v] < curLevel){
                            allAbove = false;
                            break;
                        }
                    }
                    if (allAbove)
                        ++cnt;
                }

                level_degree_map_[u] = cnt;
            }

            badVertices.push(u);

        } // end if


    } // end while

    return edge_id;
}
