#include "namespace.h"
#include "EdgeManager.h"

std::vector<EdgeIdx> EdgeManager :: getEdgeIdsAfterDuplication(std::vector<VertexIdx> &edgeV, int dupFactor)
{
    EdgeIdx startEdgeId = edgeDupList.size();
    std::vector<EdgeIdx> dupEdgeIds;
    for(int i = 0; i < dupFactor; i++)
    {
        edgeDupList.push_back(edgeV);
        dupEdgeIds.push_back(startEdgeId);
        startEdgeId += 1;
    }
    std::list<EdgeIdx>::iterator listIt = edge2IdsList[edgeV].end();
    edge2IdsList[edgeV].insert(listIt, dupEdgeIds.begin(), dupEdgeIds.end());
    return dupEdgeIds;
}

std::vector<EdgeIdx> EdgeManager :: retrieveDuplicatedEdgeIds(std::vector<VertexIdx> &edgeV, int dupFactor)
{
    // std::list<EdgeIdx> listOfEdgeIds = edge2IdsList[edgeV];
    // std::list<EdgeIdx>::iterator listIt = edge2IdsList[edgeV].begin();
    std::vector<EdgeIdx> dupEdgeIds;
    for(int i = 0; i < dupFactor; i++)
    {
        dupEdgeIds.push_back(edge2IdsList[edgeV].front());
        edge2IdsList[edgeV].pop_front();
    }
    return dupEdgeIds;
}

// int main()
// {
//     EdgeManager EM;
//     std::vector<VertexIdx> e{1,2,3,4};
//     std::vector<EdgeIdx> eidx = EM.abc(e, 1);
//     return 0;
// }