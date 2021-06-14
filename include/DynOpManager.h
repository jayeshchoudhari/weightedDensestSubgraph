#ifndef DYN_OP_MANAGER_H
#define DYN_OP_MANAGER_H

#include "namespace.h"
#include "DynamicGraph.h"
#include "EdgeManager.h"
#include <iostream>

class DynOpManager
{
    private:
        std::vector<DynamicGraph> DG;
        int active;
        int maxInstanceId;
    public:
        DynOpManager(int mId)
        {
            active= 0;
            maxInstanceId = mId;
        }
        int bindGraph(std::vector<DynamicGraph> &G);
        int addEdge(edgeVector &e, std::vector<EdgeIdx> &eDupIds, EdgeManager &EM);
};

#endif // DYN_OP_MANAGER_H