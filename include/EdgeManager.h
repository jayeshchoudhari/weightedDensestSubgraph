#ifndef EDGE_MANAGER_H
#define EDGE_MANAGER_H

#include "namespace.h"

class EdgeManager
{
    public:
        std::vector<std::vector<VertexIdx>> edgeDupList;
        std::vector<EdgeIdx> abc(std::vector<VertexIdx> edgeV, int dupFactor);
};

#endif  // EDGE_MANAGER_H