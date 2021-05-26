#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>
#include <tuple>
#include <sstream>
// #include <utility>


namespace rwNameSpace
{
    
	using VertexIdx = int64_t;
    using EdgeIdx   = int64_t;

    const VertexIdx NullVertexIdx = -1;
	const EdgeIdx NullEdgeIdx = -1;

    using Count     = int64_t;
    using Weight 	= int64_t;
    using ePair 	= std::pair<VertexIdx, VertexIdx>;
    using eTupleUnWeighted	= std::tuple<VertexIdx, VertexIdx>;
    using eTupleWeighted	= std::tuple<VertexIdx, VertexIdx, Count>;
    using VEPair 	= std::pair<VertexIdx, EdgeIdx>;
    
    using revItMapCountSetVertices = std::map<Count, std::set<VertexIdx>>::reverse_iterator;

    using edgeVector	= std::vector<VertexIdx>;
}

#endif