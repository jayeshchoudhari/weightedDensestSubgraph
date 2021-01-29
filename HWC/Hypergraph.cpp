#include "Hypergraph.hpp"
#include "Utility.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <assert.h>
#include <algorithm>

int Hypergraph::addEdge(std::vector<int>& edge) {

    std::sort(edge.begin(), edge.end());

    if (edge2ids.find(edge) != edge2ids.end() && (int)edge2ids[edge].size() >= max_keep_alive)
        return -1;

#ifdef DEBUG
    assert(edge.size() > 0);
    //make sure no duplicate vertices appear in edge
    std::unordered_set<int> _s(edge.begin(), edge.end());
    assert(_s.size() == edge.size());
#endif

    int edge_id = edgePool.size();
    edge2ids[edge].push_back(edge_id);
    edgePool.push_back(edge);
    edgePos.resize(edgePool.size());
    ++totalNumOfEdges;

    edgeCardinalityCount[(int)edge.size()]++;

    int pos = 0;
    for (const int u : edge) {
        if (edgesOf.find(u) == edgesOf.end())
            ++totalNumOfVertices;

        edgePos.back().push_back(edgesOf[u].size());
        edgesOf[u].push_back(edge_id);
        vertexPos[u].push_back(pos++);
        max_degree = std::max(max_degree, (int)edgesOf[u].size());
    }

    return edge_id;
}

int Hypergraph::removeEdge(std::vector<int>& edge) {
    std::sort(edge.begin(), edge.end());
    int edge_id = -1;
    if (edge2ids.find(edge) != edge2ids.end()){
        assert(edge2ids[edge].size() > 0);
        edge_id = edge2ids[edge].back();
        edge2ids[edge].pop_back();
        if (edge2ids[edge].empty())
            edge2ids.erase(edge);
        removeEdge(edge_id);
    }
    // this might not hold... if we allow duplicates
    // assert(edge_id != -1);
    return edge_id;
}

void Hypergraph::removeEdge(const int edge_id) {

#ifdef DEBUG
    assert (edge_id < (int)edgePool.size());
    assert (edgePool[edge_id].size() > 0);
#endif

    int cardinality = (int)edgePool[edge_id].size();
    assert(cardinality >= 1);
    edgeCardinalityCount[cardinality]--;
    if (edgeCardinalityCount[cardinality] == 0)
        edgeCardinalityCount.erase(cardinality);

    for (int i = 0; i < (int)edgePool[edge_id].size(); i++) {
        int u = edgePool[edge_id][i];
        int pos = edgePos[edge_id][i];

        assert(edgesOf[u].size());
        if (pos != (int)edgesOf[u].size() - 1) {
            int last_edge_id = edgesOf[u].back();
            int last_pos = vertexPos[u].back();
            edgePos[last_edge_id][last_pos] = pos;

            std::swap(edgesOf[u][pos], edgesOf[u].back());
            std::swap(vertexPos[u][pos], vertexPos[u].back());
        }

        edgesOf[u].pop_back();
        vertexPos[u].pop_back();

// By default
// do not remove vertex even when it's degree is 0
if(REMOVE_VERTEX_WHEN_EDGE_EMPTY){
        if (edgesOf[u].size() == 0){
            edgesOf.erase(u);
            vertexPos.erase(u);
            --totalNumOfVertices;
        }
}
    
    }
        // release memory
    edgePool[edge_id].clear();
    edgePos[edge_id].clear();
    --totalNumOfEdges;
        
    return;
}

void Hypergraph::removeVertex(const int u) {
    while (edgesOf[u].size() > 0) {
        removeEdge(edgesOf[u].back());    
    }
    edgesOf.erase(u);
    vertexPos.erase(u);
    --totalNumOfVertices;
}



int Hypergraph::degree(const int u) {
    if (edgesOf.find(u) == edgesOf.end())
        return 0;
    
    return (int)(edgesOf[u].size());
}

int Hypergraph::numOfVertices() const {
    assert(totalNumOfVertices == (int)edgesOf.size());

    return totalNumOfVertices;
}

int Hypergraph::numOfEdges() const {
    return totalNumOfEdges;
}

// int Hypergraph::degreeUpperBound() const {
//     return max_cardinality;
// }

bool Hypergraph::contains(int u)  {
    return edgesOf.find(u) != edgesOf.end();
}

void Hypergraph::fillDegree(std::vector<std::vector<int>>& bin, std::unordered_map<int, int>& pos) {
    for (auto it = edgesOf.begin(); it != edgesOf.end(); it++) {
        int u = it->first;
        int d = degree(u);
        while((int)bin.size() < d+1)
            bin.push_back(std::vector<int>());
        pos[u] = bin[d].size();
        bin[d].push_back(u);
    }
}

void Hypergraph::neighborsOf(const int u, std::vector<int> * vec){
    std::unordered_set<int> s;
    for (int edge_id : edgesOf[u])
        for (int v : edgePool[edge_id])
            s.insert(v);
    for(int v : s)
        vec->push_back(v);
}

void Hypergraph::edgesContain(const int u, std::vector<int>* vs){
    if (this->contains(u))
        for (int v : edgesOf[u])
            vs->push_back(v);
    return;
}

void Hypergraph::edgeVertices(int edge_id, std::vector<int> * vp){
    assert(edge_id < (int)edgePool.size());
    for (const int v :  edgePool[edge_id])
        vp->push_back(v);
    return;
}

void Hypergraph::vertices(std::vector<int> * vp) const {
    for (auto iter = edgesOf.begin(); iter != edgesOf.end(); ++iter)
        vp->push_back(iter->first);
}


double Hypergraph::graphDensity() const {
    if (totalNumOfVertices == 0)
        return 0.0;
    return (1.0*totalNumOfEdges)/totalNumOfVertices;
}
void Hypergraph::show(){
    printVecVec(edgePool);
}

int Hypergraph::maxEdgeCardinality() const {
    if (edgeCardinalityCount.empty())
        return 2;
    // assert (edgeCardinalityCount.rbegin()->first >= 2);
    if (max_cardinality_threshold != -1) {
        return min(max_cardinality_threshold, (int)edgeCardinalityCount.rbegin()->first);
    }
    return edgeCardinalityCount.rbegin()->first;
}

void Hypergraph::edges(vector<int> * vec) {
    unordered_set<int> edges;

    for (auto it = edgesOf.begin(); it != edgesOf.end(); it++) {
        for (const int edge_id : it->second) {
            edges.insert(edge_id);
        }
    }

    for (const int edge_id : edges) {
        vec->push_back(edge_id);
    }

    return;
}