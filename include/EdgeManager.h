#ifndef EDGE_MANAGER_H
#define EDGE_MANAGER_H

#include "namespace.h"

class EdgeManager
{
    public:
        EdgeIdx currentEdgeId;
        std::vector<std::vector<VertexIdx>> edgeDupList;
        std::unordered_map<EdgeIdx, std::vector<VertexIdx>> edgeDupMap;
        std::unordered_map<std::vector<VertexIdx>, std::list<EdgeIdx>, vectorHash> edge2IdsList;

        EdgeManager()
        {
            currentEdgeId = 0;
        }
        std::vector<EdgeIdx> getEdgeIdsAfterDuplication(std::vector<VertexIdx> &edgeV, Count dupFactor);
        std::vector<EdgeIdx> retrieveDuplicatedEdgeIds(std::vector<VertexIdx> &edgeV, Count dupFactor);
        int removeEdgeFromMemory(std::vector<EdgeIdx> edgeDuplicatorIds);
};

#endif  // EDGE_MANAGER_H