#ifndef NAMESPACE_H
#define NAMESPACE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <queue>
#include <algorithm>
#include <string>
#include <set>
#include <list>

using namespace std;

using VertexIdx = int64_t;
using EdgeIdx   = int64_t;
using Count     = int64_t;
using Weight 	= int64_t;
using ePair 	= std::pair<VertexIdx, VertexIdx>;
using edgeVector	= std::vector<VertexIdx>;

using VEPair 	= std::pair<VertexIdx, EdgeIdx>;
using eTupleUnWeighted	= std::tuple<VertexIdx, VertexIdx>;
using eTupleWeighted	= std::tuple<VertexIdx, VertexIdx, Count>;

using revItMapCountSetVertices = std::map<Count, std::set<VertexIdx>>::reverse_iterator;

const VertexIdx NullVertexIdx = -1;
const EdgeIdx NullEdgeIdx = -1;


#endif      // NAMESPACE_H