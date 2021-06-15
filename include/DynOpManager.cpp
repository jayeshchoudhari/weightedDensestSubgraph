#include "namespace.h"
#include "DynOpManager.h"
#include "EdgeManager.h"
#include "DynamicGraph.h"
#include <iostream>

// using namespace std;

DynOpManager :: DynOpManager(int maxId)
{
    active= 0;
    maxInstanceId = maxId;
}

int DynOpManager :: bindGraph(std::vector<DynamicGraph> &G)
{
    DG = G;
    return 0;
}

int DynOpManager :: addEdge(edgeVector &currentEdge, std::vector<EdgeIdx> &eDupId, EdgeManager &EM)
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


int DynOpManager :: removeEdge(edgeVector &currentEdge, std::vector<EdgeIdx> &eDupId, EdgeManager &EM)
{
    double rho, rhoEst, rhoEstActive;
    for(int dup = 0; dup < eDupId.size(); ++dup)
    {
        // EdgeIdx delEdgeId = DG[0].checkEdgeExistence(currentEdge);
        EdgeIdx delEdgeId = eDupId[dup];
        if(delEdgeId != NullEdgeIdx)
        {
            // edgeDeletions += 1;
            // totalProcessedEdges += 1;

            // beginClock = std::chrono::steady_clock::now();
            for(int copy = maxInstanceId; copy >= active; --copy)
            {
                int checkEdge = DG[copy].checkEdgeExistenceInPendingList(delEdgeId);
                if(checkEdge != 0)
                {
                    DG[copy].removeEdgeFromPendingList(currentEdge, delEdgeId);
                }
                else
                {
                    VertexIdx lastVertex = DG[copy].deleteEdge(currentEdge, delEdgeId, EM);
                }
            }

            if(active - 1 >= 1)
            {
                rho = DG[active].getDensity();
                rhoEstActive = DG[active].getRhoEst();

                if(rho < rhoEstActive)
                {
                    active = active - 1;
                    VertexIdx lastVertex = DG[active].deleteEdge(currentEdge, delEdgeId, EM);
                }
            }

            for(int copy = active - 1; copy >= 1; --copy)
            {
                int checkEdge = DG[copy].checkEdgeExistenceInPendingList(currentEdge, delEdgeId);

                if(checkEdge != 0)
                {
                    DG[copy].removeEdgeFromPendingList(currentEdge, delEdgeId);
                }
                else
                {
                    // VertexIdx lastVertex = DG[copy].deleteEdgeReturnLastVertex(currentEdge, delEdgeId);
                    VertexIdx lastVertex = DG[copy].deleteEdge(currentEdge, delEdgeId, EM);
                    // std::pair<EdgeIdx, edgeVector> eIdEdgePair = DG[copy].getPendingEdgeForLastVertex(lastVertex);
                    EdgeIdx eId = DG[copy].getPendingEdgeForLastVertex(lastVertex);
                    if(eId != NullEdgeIdx)
                    {
                        edgeVector e = EM.edgeDupList[eId];
                        DG[copy].insertEdge(e, eId, EM);
                        DG[copy].removeEdgeFromPendingList(e, eId);
                    }
                }
            }

            // DG[0].removeEdgeFromMap(currentEdge);

            endClock = std::chrono::steady_clock::now();
            localDeleteTime = std::chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
            perReportTime += localDeleteTime;
            /*
            */
        }
        else
        {
            std::cout << "Edge does not exists to delete...\n";
            break;
        }
    }
    return 0;
}