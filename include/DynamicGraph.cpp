#include "namespace.h"
#include "DynamicGraph.h"
#include "GraphLoader.h"
#include <iostream>

using namespace std;

int DynamicGraph :: initializeVariables(int i, Count nv, float epsUD)
{
    nVertices = nv;
	epsVal = epsUD;
    alpha = round(1.0/(pow(epsVal,2)));

	if(i > 0)
    {
        rhoEst = pow(2, i-2) * alpha;
        eta = (2 * rhoEst)/alpha;
    }
    else
    {
        rhoEst = 0;
        eta = 0;
    }
	
    // // Maintain list of InNbrs...
    InNbrs.resize(nv);

    // // ********* TO MAINTAIN A UPDATING PRIORITY QUEUE OF OUTNEIGHBORS *************** 
    // // this is to be done for each node.. and thus a vector here...
    // // Each element of a vector is a map --  where the map is with a key Count, and value is a set of 
    // // of edgeIds, which have the headVertices with value (indegree) as Count...
    // // and nodeToOutdegMap is a reverseMap of the same...
	outdegToNodeMap.resize(nv);
	nodeToOutdegMap.resize(nv);
	InDegreeFromNodesView.resize(nv);


	
    // // to maintain visitNext structure... 
	listOfNeighbors.resize(nv);
	mapToNeighborsList.resize(nv);
	nextPositionIteratorInc.resize(nv);
	nextPositionIteratorDec.resize(nv);
	nextPositionIteratorTightInNbr.resize(nv);

	return 0;
}
