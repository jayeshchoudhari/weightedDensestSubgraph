#include "namespace.h"
#include "EdgeManager.h"

std::vector<EdgeIdx> EdgeManager :: abc(std::vector<VertexIdx> edgeV, int dupFactor)
{
    EdgeIdx startEdgeId = edgeDupList.size();
    std::vector<EdgeIdx> dupEdgeIds;
    for(int i = 0; i < dupFactor; i++)
    {
        edgeDupList.push_back(edgeV);
        dupEdgeIds.push_back(startEdgeId);
        startEdgeId += 1;
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