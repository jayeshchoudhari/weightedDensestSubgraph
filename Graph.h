

#ifndef rwNameSpace
#define rwNameSpace

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

using namespace std;

namespace rwNameSpace {

	using VertexIdx = int64_t;
    using EdgeIdx   = int64_t;

    const VertexIdx NullVertexIdx = -1;
	const EdgeIdx NullEdgeIdx = -1;

    using Count     = int64_t;
    using Weight 	= int64_t;
    using ePair 	= pair<VertexIdx, VertexIdx>;
    using eTupleUnWeighted	= tuple<VertexIdx, VertexIdx>;
    using eTupleWeighted	= tuple<VertexIdx, VertexIdx, Count>;
    using VEPair 	= pair<VertexIdx, EdgeIdx>;
    
    using revItMapCountSetVertices = map<Count, set<VertexIdx>>::reverse_iterator;

    using edgeVector	= vector<VertexIdx>;
}

#endif