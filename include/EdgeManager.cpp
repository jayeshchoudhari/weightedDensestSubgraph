#include "namespace.h"
#include "EdgeManager.h"

std::vector<EdgeIdx> EdgeManager :: getEdgeIdsAfterDuplication(std::vector<VertexIdx> &edgeV, Count dupFactor)
{
    std::vector<EdgeIdx> dupEdgeIds;
    for(Count i = 0; i < dupFactor; i++)
    {
        // edgeDupList.push_back(edgeV);
        dupEdgeIds.push_back(currentEdgeId);
        edgeDupMap[currentEdgeId] = edgeV;
        currentEdgeId += 1;
    }
    std::list<EdgeIdx>::iterator listIt = edge2IdsList[edgeV].end();
    edge2IdsList[edgeV].insert(listIt, dupEdgeIds.begin(), dupEdgeIds.end());
    return dupEdgeIds;
}

std::vector<EdgeIdx> EdgeManager :: retrieveDuplicatedEdgeIds(std::vector<VertexIdx> &edgeV, Count dupFactor)
{
    std::vector<EdgeIdx> dupEdgeIds;
    for(Count i = 0; i < dupFactor; i++)
    {
        dupEdgeIds.push_back(edge2IdsList[edgeV].front());
        edge2IdsList[edgeV].pop_front();
    }
    return dupEdgeIds;
}

int EdgeManager :: removeEdgeFromMemory(std::vector<EdgeIdx> edgeDuplicatorIds)
{
    for(Count i = 0; i < edgeDuplicatorIds.size(); i++)
    {
        edgeDupMap.erase(edgeDuplicatorIds[i]);
    }

    return 0;
}

// int main()
// {
//     EdgeManager EM;
//     std::vector<VertexIdx> e{1,2,3,4};
//     std::vector<EdgeIdx> eidx = EM.abc(e, 1);
//     return 0;
// }