#ifndef EDGE_MANAGER_H
#define EDGE_MANAGER_H

#include "namespace.h"

class EdgeManager
{
    public:
        std::vector<std::vector<VertexIdx>> edgeDupList;
        std::unordered_map<std::vector<VertexIdx>, std::list<EdgeIdx>, vectorHash> edge2IdsList;
        std::vector<EdgeIdx> getEdgeIdsAfterDuplication(std::vector<VertexIdx> &edgeV, int dupFactor);
        std::vector<EdgeIdx> retrieveDuplicatedEdgeIds(std::vector<VertexIdx> &edgeV, int dupFactor);
};

#endif  // EDGE_MANAGER_H