#include "namespace.h"
#include "DynOpManager.h"
#include "EdgeManager.h"
#include "DynamicGraph.h"
#include <iostream>

// using namespace std;

DynOpManager :: DynOpManager(int maxId, std::string file_name)
{
    active= 0;
    maxInstanceId = maxId;
    outFile.open(file_name.c_str(), std::ios::out);
    outFile << "Label  MaxInDeg  OneMinusEpsMaxInDeg  MaxDensity  Size  TimePerWindow(micros)  NumOpPerWindow  TimePerOp\n";
    outFile.flush();
}

int DynOpManager :: bindGraph(std::vector<DynamicGraph> &G)
{
    DG = G;
    return 0;
}

int DynOpManager :: addEdge(edgeVector &currentEdge, std::vector<EdgeIdx> &eDupId, EdgeManager &EM)
{
    double rho, rhoEst, rhoEstActive;
    // Duplication of edges...
    for(unsigned int dup = 0; dup < eDupId.size(); ++dup)
    {
        // Affordable copies...
        for(int copy = maxInstanceId; copy >= active + 1; --copy)
        {
            DG[copy].insertEdge(currentEdge, eDupId[dup], EM);
        }

        if(active + 1 <= maxInstanceId)
        {
            rho = DG[active + 1].getDensity();
            rhoEstActive = DG[active].getRhoEst();      // change this -- not to call the function -- make rhoEst public...

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
            // std::cout  << "updating all instances --- \n";
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
            // std::cout  << "setting up active instance --- \n";
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
            // std::cout  << "checking for pending edges...--- \n";
            for(int copy = active - 1; copy >= 1; --copy)
            {
                int checkEdge = DG[copy].checkEdgeExistenceInPendingList(delEdgeId);

                if(checkEdge != 0)
                {
                    // std::cout << "Exists and Removing from pending list --- \n";
                    DG[copy].removeEdgeFromPendingList(currentEdge, delEdgeId);
                }
                else
                {
                    // std::cout << copy << " -- Not in pending list -- Getting the last vertex after deletion --- \n";
                    // VertexIdx lastVertex = DG[copy].deleteEdgeReturnLastVertex(currentEdge, delEdgeId);
                    VertexIdx lastVertex = DG[copy].deleteEdge(currentEdge, delEdgeId, EM);
                    // std::pair<EdgeIdx, edgeVector> eIdEdgePair = DG[copy].getPendingEdgeForLastVertex(lastVertex);
                    
                    /*
                    std::cout << "Getting Pending edges of the last vertex--- \n";
                    EdgeIdx eId = DG[copy].getPendingEdgeForLastVertex(lastVertex);
                    if(eId != NullEdgeIdx)
                    {
                        edgeVector e = EM.edgeDupMap[eId];
                        std::cout << "Adding edge to the instance --- \n";

                        for(int i = 0; i < e.size(); i++)
                            std::cout << e[i] << " ";
                        std::cout << "\n";
                        
                        DG[copy].insertEdge(e, eId, EM);
                        std::cout << "Removing from pending list.. --- \n";
                        DG[copy].removeEdgeFromPendingList(e, eId);
                    }
                    */
                    // std::cout << "Getting Pending edges of the last vertex--- \n";
                    std::pair<EdgeIdx, edgeVector> eIdEdgePair = DG[copy].getPendingEdgeForLastVertex(lastVertex);

                    if(eIdEdgePair.first != NullEdgeIdx)
                    {
                        // std::cout << eIdEdgePair.first << " --- Adding edge to the instance --- \n";
                        // for(int i = 0; i < eIdEdgePair.second.size(); i++)
                        //     std::cout << eIdEdgePair.second[i] << " ";
                        // std::cout << "\n";
                        DG[copy].insertEdge(eIdEdgePair.second, eIdEdgePair.first, EM);
                        // std::cout << "Removing from pending list.. --- \n";
                        DG[copy].removeEdgeFromPendingList(eIdEdgePair.second, eIdEdgePair.first);
                    }

                }
            }

            // DG[0].removeEdgeFromMap(currentEdge);

            // endClock = std::chrono::steady_clock::now();
            // localDeleteTime = std::chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
            // perReportTime += localDeleteTime;
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

int DynOpManager :: addPendingEdgesToActiveInstance(EdgeManager &EM)
{
    // if(DG[DOM.active].getPendingCount() > 0)
    if(DG[active].pendingListOfEdges.size() > 0)
    {
        // std::cout << "I am coming here.." << DG[active].getPendingCount() << "\n";
        DG[active].insertListOfPendingEdges(EM);
        // std::cout << "I am coming here.." << DG[active].getPendingCount() << "\n";
    }
    return 0;
}

int DynOpManager :: getDensityEstimate(EdgeManager &EM, double localAlpha, vectorListMap &mainEdge2Ids, double micros, Count numOpPerWindow, std::string label)
{
    // Count yearLabel;
    // ss >> yearLabel;

    // std::cout << "Total Processed Edges -- " << totalProcessedEdges << "\n";
    // std::cout << "Edge Additions = " << edgeAdditions << " Edge Deletions = " << edgeDeletions << "\n";

    /*
    std::vector<double> OneMinusEpsMaxInDegVector;
    std::vector<double> MaxInDegVector;

    MaxInDegVector.push_back(0);
    OneMinusEpsMaxInDegVector.push_back(0);

    for(int copy = 1; copy <= maxInstanceId; ++copy)
    {
        maxInDeg = (DG[copy].getMaxLabel()*1.0)/localAlpha;
        MaxInDegVector.push_back(maxInDeg);
        OneMinusEpsMaxInDeg = DG[copy].getDensity()/localAlpha;
        OneMinusEpsMaxInDegVector.push_back(OneMinusEpsMaxInDeg);
    }
    */
    double maxInDeg, OneMinusEpsMaxInDeg;
    maxInDeg = (DG[active].getMaxLabel()*1.0)/localAlpha;
    OneMinusEpsMaxInDeg = DG[active].getDensity()/localAlpha;

    Count duplicationFactor = (Count)localAlpha;
    std::set<VertexIdx> denseSubgraph = DG[active].querySubgraph(localAlpha);
    unsigned int denseSubgraphSize = denseSubgraph.size();
    double estimatedDensity =  DG[active].getDensityOfInducedSubgraph(denseSubgraph, mainEdge2Ids, duplicationFactor)/localAlpha;

    std::cout << "Getting max partition density....\n";
    std::pair<double, unsigned int> maxPartitionedDensity = DG[active].getMaxPartitionDensity(mainEdge2Ids, duplicationFactor);

    std::cout << "Reporting Density... - " << maxPartitionedDensity.first << " " << maxPartitionedDensity.second << std::endl;

    outFile << label << " " << maxInDeg << " " << OneMinusEpsMaxInDeg << " " << maxPartitionedDensity.first << " " << maxPartitionedDensity.second << " " << micros << " " << numOpPerWindow << " " << micros/numOpPerWindow << "\n";
    
    return 0;
}