#ifndef DYN_OP_MANAGER_H
#define DYN_OP_MANAGER_H

#include "namespace.h"
#include "DynamicGraph.h"
#include "EdgeManager.h"
#include <iostream>

class DynOpManager
{
    private:
        int active;
        int maxInstanceId;
        std::vector<DynamicGraph> DG;
    public:
        DynOpManager(int maxId);
        int bindGraph(std::vector<DynamicGraph> &G);
        int addEdge(edgeVector &e, std::vector<EdgeIdx> &eDupIds, EdgeManager &EM);
        int removeEdge(edgeVector &currentEdge, std::vector<EdgeIdx> &eDupId, EdgeManager &EM);
        int getDensityEstimate(EdgeManager &EM, double localAlpha, vectorListMap &mainEdge2Ids);
};

#endif // DYN_OP_MANAGER_H