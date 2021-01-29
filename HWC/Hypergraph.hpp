#ifndef __HYPERGRAPH__
#define __HYPERGRAPH__

#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>

struct vectorHash{
    std::size_t operator()(std::vector<int> const& vec) const {
      std::size_t seed = vec.size();
      for(auto& i : vec) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
};


class Hypergraph{


public:
    int numOfVertices() const;
    int numOfEdges() const;

    int addEdge(std::vector<int>& edge);
    int removeEdge(std::vector<int>& edge);
    void removeVertex(const int vertex_id);

    void neighborsOf(const int u, std::vector<int> * vec);
    void edgesContain(const int u, std::vector<int>* vs);
 

    int degree(const int u);
    // int degreeUpperBound() const;
    int maxEdgeCardinality() const;
    int maxDegree() const {
        if (max_degree_threshold != -1) {
            return std::min(max_degree, max_degree_threshold);
        }
        return max_degree;    
    }

    void fillDegree(std::vector<std::vector<int>>& bin, std::unordered_map<int, int>& pos);
    void vertices(std::vector<int> * vp) const;

    bool contains(int u) ;
    double graphDensity() const;
    void edgeVertices(int edge_id, std::vector<int> * vec);
    void edges(std::vector<int> * vec);

    void show();

    Hypergraph(bool rv, int edge_card = -1, int m_degree = -1) {
        totalNumOfVertices = totalNumOfEdges = 0;
        max_degree = max_cardinality = 0;
        REMOVE_VERTEX_WHEN_EDGE_EMPTY = rv;
        max_cardinality_threshold = edge_card;
        max_degree_threshold = m_degree;
    }

    ~Hypergraph(){
        // std::cerr << " Hypergraph :" << std::endl 
        //           << " max_cardinality: " << max_cardinality 
        //           << "max degree :" << max_degree << std::endl;
    }

private:
    void removeEdge(const int edge_id);
    
    bool REMOVE_VERTEX_WHEN_EDGE_EMPTY;

    static const int max_keep_alive = 1;

    int totalNumOfVertices;
    int totalNumOfEdges;


    std::map<int, int> edgeCardinalityCount;

    int max_cardinality;
    int max_degree;

    int max_cardinality_threshold;
    int max_degree_threshold;


    std::vector<std::vector<int>> edgePool;
    std::vector<std::vector<int>> edgePos;

    std::unordered_map<int, std::vector<int>> edgesOf;
    std::unordered_map<int, std::vector<int>> vertexPos;

    // each hyperege can appear at most max_keep_alive times
    std::unordered_map<std::vector<int>, std::vector<int>, vectorHash> edge2ids;
};

#endif