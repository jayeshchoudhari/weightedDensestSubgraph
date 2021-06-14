#include "namespace.h"
#include "DynOpManager.h"
#include "EdgeManager.h"
#include "DynamicGraph.h"
#include <iostream>

using namespace std;

int DynOpManager :: bindGraph(vector<DynamicGraph> &G)
{
    DG = G;
    return 0;
}

int DynOpManager :: addEdge(edgeVector &currentEdge, vector<EdgeIdx> &eDupId, EdgeManager &EM)
{
    double rho, rhoEst, rhoEstActive;
    for(unsigned int dup = 0; dup < eDupId.size(); ++dup)
    {

        for(int copy = maxInstanceId; copy >= active + 1; --copy)
        {
            DG[copy].insertEdge(currentEdge, eDupId[dup], EM);
        }

        if(active + 1 <= maxInstanceId)
        {
            rho = DG[active + 1].getDensity();
            rhoEstActive = DG[active].getRhoEst();

            if(rho >= 2*rhoEstActive)
            {
                active = active + 1;
            }
            else
            {
                DG[active].insertEdge(currentEdge, eDupId[dup], EM);
            }
        }
        else
        {
            DG[active].insertEdge(currentEdge, eDupId[dup], EM);
        }

        for(int copy = active - 1; copy >= 1; --copy)
        {
            rhoEst = DG[copy].getRhoEst();
            Count minLoad = DG[copy].getMinLoadInE(eDupId[dup], EM);

            if(minLoad >= 2 * rhoEst)
            {
                DG[copy].addEdgeToPendingList(currentEdge, eDupId[dup]);
            }
            else
            {
                DG[copy].insertEdge(currentEdge, eDupId[dup], EM);
            }
        }
    }
    return 0;
}


