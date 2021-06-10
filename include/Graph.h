

#ifndef rwNameSpace
#define rwNameSpace

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <fstream>

#include <vector>
#include <unordered_map>
#include <queue>
// #include <utility>

using namespace std;

namespace rwNameSpace {

	using VertexIdx = int64_t;
    using EdgeIdx   = int64_t;
    using Count     = int64_t;
    using Weight 	= int64_t;
    using ePair 	= std::pair<VertexIdx, VertexIdx>;
    const EdgeIdx invalidEdge = -1;
}

#endif