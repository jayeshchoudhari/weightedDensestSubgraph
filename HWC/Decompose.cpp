#include "ApproximateDSFullyDyn.hpp"
#include "Decompose.hpp"
#include "Utility.hpp"
#include "Hypergraph.hpp"

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


void Decompose::updateTao(){
    assert(alpha_beta_ >= 1);
    tao_ = max_iter((graph_ == NULL ? 0 : graph_->numOfVertices()), epsilon_);
}

void Decompose::Check(vector<int>& edge, int le){    

    stack<int> bad;

    for (const int u : edge) {
        if (level_map_[u] == le){
            bad.push(u);
        }
    }

    updateTao();
    touchedLevel.clear();
    touchedLevel.insert(le);
    
    while (!bad.empty()){
        int u = bad.top();
        bad.pop();
        int lu = level_map_[u];
        int curDegree = (int)rmEdge_map_[make_pair(u, lu)].size();
        int preDegree;
        if (lu > 1){
            preDegree = curDegree + (int)rmEdge_map_[make_pair(u, lu-1)].size();
        }

        if (lu < tao_ && curDegree >= alpha_beta_) {
            Promote(u, bad);
        } else if (lu > 1 && preDegree < beta_){
            Demote(u, bad);
        } 
    }

    updateDensity();
}

void Decompose::Add(int u, int l, int edge_id){
    pair<int, int> ul = make_pair(u, l);
    rmEdge_map_[ul].push_back(edge_id);
    return;
}

void Decompose::Erase(int u, int l, int edge_id){
    pair<int, int> pii = make_pair(u, l);
    vector<int>::iterator it = find(rmEdge_map_[pii].begin(), rmEdge_map_[pii].end(), edge_id);
//    std::cout <<"u=" << u <<" l=" << l << " Edge id = " << edge_id << std::endl;
    if (it == rmEdge_map_[pii].end())
    {
        std::cout <<"Erase:"<<"u=" << u <<" l=" << l << " Edge id = " << edge_id << std::endl;
    }
    assert(it != rmEdge_map_[pii].end());
    swap(*it, rmEdge_map_[pii].back());
    rmEdge_map_[pii].pop_back();
    return;
}

void Decompose::Promote(int u, stack<int>& bad) {
    int t = level_map_[u];
    level_map_[u] = t + 1;
    level_vertices_[t+1] += 1;
    pair<int, int> ut = make_pair(u, t);

    touchedLevel.insert(t);
    touchedLevel.insert(t+1);

    vector<int> & vi = rmEdge_map_[ut];
    int i = 0; 
    while (i < (int)vi.size()) {
        int edge_id = vi[i];
        int min_level = INF_LEVEL;
        vector<int> vertices;
        graph_-> edgeVertices(edge_id, &vertices);
        assert(vertices.size() >= 1);
        bool earlyBreak = false;
        for (const int v : vertices) {
            assert(level_map_[v] >= t);
            min_level = min(min_level, level_map_[v]);
            if (min_level < t + 1){
                earlyBreak = true;
                break;
            }
        }

        if (earlyBreak){
            assert(min_level == t);
            i++;
            continue;
        }

        assert(min_level == t+1);
        for (const int v : vertices) {
            if (u != v) {
                Erase(v, t, edge_id);
                Add(v, t+1, edge_id);
                
                if (level_map_[v] == t+1){
                    bad.push(v);
                }

            } else {
                Add(u, t+1, edge_id);
                swap(vi[i], vi.back());
                vi.pop_back();
            }
        }

        level_edges_[t+1] += 1;
    }

    bad.push(u);

}


void Decompose::Demote(int u, stack<int>& bad) {
    int t = level_map_[u];
    pair<int, int> ut = make_pair(u, t);
    pair<int, int> ut1 = make_pair(u, t-1);
    
    touchedLevel.insert(t);
    touchedLevel.insert(t-1);
    
    for (int edge_id : rmEdge_map_[ut]) {
        vector<int> vertices;
        graph_-> edgeVertices(edge_id, &vertices);
        int min_level = INF_LEVEL;
        for (const int v : vertices) {
            if (v != u){
                assert(level_map_[v] >= t);
                min_level = min(min_level, level_map_[v]);
            }
        }
        assert(min_level >= t);
        for (const int v : vertices) {
            if (u != v){
                Erase(v, t, edge_id);
                Add(v, t-1, edge_id);
                if (level_map_[v] == t){
                    bad.push(v);
                }
            }
        }
        --level_edges_[t];            
    }


    level_map_[u] = t - 1;
    level_vertices_[t] -= 1;

    for (int edge_id : rmEdge_map_[ut])
        rmEdge_map_[ut1].push_back(edge_id);
    rmEdge_map_.erase(ut);

    bad.push(u);
}


double Decompose::addEdge(vector<int>& edge, int edge_id) {

    // cout << "adding " << beta_ << " | " << alpha_beta_ << " : " << edge_id << endl;
    // printLnVec(edge);
    
    int min_level = INF_LEVEL;
    for (const int u : edge) {
        if (level_map_.find(u) == level_map_.end()) {
            level_map_[u] = 1;
            level_vertices_[1] += 1;
        }
        min_level = min(min_level, level_map_[u]);
    }
    
    assert(min_level != INF_LEVEL);

    for (const int u : edge){
        Add(u, min_level, edge_id);
    }

    for (int i = 1; i <= min_level; i++)
        level_edges_[i] += 1;

    Check(edge, min_level);

    return density();
}

double Decompose::removeEdge(vector<int>& edge, int edge_id) { 
    int min_level = INF_LEVEL;
    for (const int u : edge) {
        assert(level_map_.find(u) != level_map_.end());
        min_level = min(min_level, level_map_[u]);
    }
    assert(min_level != INF_LEVEL);

    for (const int u : edge) {
        Erase(u, min_level, edge_id);
    }

    for (int i = 1; i <= min_level; i++)
        level_edges_[i] -= 1;
    
    Check(edge, min_level);

    return density();
}

void Decompose::newConstruct(){
    vector<int> edges;
    graph_->edges(&edges);
    for (int edge_id : edges){
        vector<int> vertices;
        graph_->edgeVertices(edge_id, &vertices);
        addEdge(vertices, edge_id);
    }
    return;
}

void Decompose::updateDensity(){
    for (int l : touchedLevel){
        double density;
        if(level_vertices_[l] == 0) {
            density = 0.0;
        } else {
            density = (1.0 * level_edges_[l]) / level_vertices_[l];
        } 
        
        if (level_density_.find(l) == level_density_.end()) {
            if (density > 0.0)    {
                level_density_set_.insert(make_pair(density, l));
                level_density_[l] = density;
            }
        } else {
            double oldDensity = level_density_[l];
            if (oldDensity != density){
                int tmp = level_density_set_.erase(make_pair(oldDensity, l));
                assert(tmp == 1);
                level_density_set_.insert(make_pair(density, l));
                level_density_[l] = density;
            }
        }
    }
}

int Decompose::subgraphSize(){
    if (level_vertices_.empty()) {
        return 0;
    }
    int bestLevel = (*level_density_set_.rbegin()).second;
    return level_vertices_[bestLevel];
}

// pair<int, int> Decompose::subgraph(){
//     if (level_density_set_.empty())
//         return make_pair(0, 0);
//     int bestLevel = (*level_density_set_.rbegin()).second;
//     return make_pair(level_vertices_[bestLevel], level_edges_[bestLevel]);
// }

double Decompose::density(){
    if (level_density_set_.empty())
        return 0.0;
    int bestLevel = (*level_density_set_.rbegin()).second;
    return (1.0 * level_edges_[bestLevel])/level_vertices_[bestLevel];
}





