#ifndef DYN_OP_MANAGER_H
#define DYN_OP_MANAGER_H

#include "namespace.h"
#include "DynamicGraph.h"
#include "EdgeManager.h"
#include <iostream>

class DynOpManager
{
    private:
        int maxInstanceId;
        std::vector<DynamicGraph> DG;
        std::ofstream outFile;
    public:
        int active;
        DynOpManager(int maxId, std::string file_name);
        int bindGraph(std::vector<DynamicGraph> &G);
        int addEdge(edgeVector &e, std::vector<EdgeIdx> &eDupIds, EdgeManager &EM);
        int removeEdge(edgeVector &currentEdge, std::vector<EdgeIdx> &eDupId, EdgeManager &EM);
        int addPendingEdgesToActiveInstance(EdgeManager &EM);
        int getDensityEstimate(EdgeManager &EM, double localAlpha, vectorListMap &mainEdge2Ids, double micros, Count numOpPerWindow, std::string label);
};

#endif // DYN_OP_MANAGER_H