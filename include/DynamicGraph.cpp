#include "namespace.h"
#include "DynamicGraph.h"
#include "GraphLoader.h"
#include <iostream>

using namespace std;

DynamicGraph :: DynamicGraph(int i, Count nv, float epsUD)
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
}

int DynamicGraph :: insertEdge(edgeVector &currentEdge, EdgeIdx eId, EdgeManager &EM)
{
	// // e = u,v -- at this moment an edge isnt directed...
	// // hereafter it will be... from here onwards we would orient the edges...

	VertexIdx w, wPrime;
	EdgeIdx ePrime, lastEId;

	w = getMinDegreeVertexInE(currentEdge);
	lastEId = eId;

	// // std::cout << nodeInDeg[w] << " headNode indegree during addition\n";
	addDirectedEdgeToInOutNbrs(currentEdge, eId, w);

	// // std::cout << InNbrs[w].size() << " headNode innbrs size during addition\n";

	// // headOfEdgeId[eId] = w;

	// // check if this results into making some neighboring edge of w tight...
	std::pair<VertexIdx, EdgeIdx> minNbrNodeEdgePair = getTightInNbr(w, EM);
	ePrime = minNbrNodeEdgePair.second;		// tight edge... -- head of this edge give tight In-neighbor node...
	// // i.e. an edge with a headNode whose degree is very less as compared to that of w....
	while(ePrime != -1)
	{
		wPrime = minNbrNodeEdgePair.first;
		// // wPrime's in-deg is very less...
		// // so flip the edge to wPrime..?

		// // eTupleUnWeighted eFlip(wPrime, w);
		flipDirectedEdge(ePrime, w, wPrime, EM);		// flipDirectedEdge(eId, oldHeadNode, newHeadNode)
		// // so now wPrime is settled -- i.e. we have increased its indegree 
		// // now we need to check if any in-neighbor of wPrime, violates the condition or has very less indegree...
		w = wPrime;
		lastEId = ePrime;
		minNbrNodeEdgePair = getTightInNbr(w, EM);
		ePrime = minNbrNodeEdgePair.second;
	}

	incrementDu(w, EM);

	return 0;
}


VertexIdx DynamicGraph :: deleteEdge(edgeVector &currentEdge, EdgeIdx eId, EdgeManager &EM)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	VertexIdx w, wPrime;
	EdgeIdx ePrime, lastEId;

	VertexIdx headNode;

	// if(headOfEdgeId.find(eId) != headOfEdgeId.end())
	// {
	headNode = headOfEdgeId[eId];
	// }
	// else
	// {
	// 	headNode = NullVertexIdx;
	// }

	// std::cout << nodeInDeg[headNode] << " headNode indegree before deletion\n";
	// std::cout << InNbrs[headNode].size() << " headNode innbrs size before deletion\n";
	
	// if(headNode != NullVertexIdx)
	// {
		// remove e from the InNbrs of headNode...
	removeDirectedEdgeFromInOutNbrs(currentEdge, eId, headNode);

	lastEId = eId;
	w = headNode;
	ePrime = getTightOutNbr(w);
	
	std::map<EdgeIdx, int> flippedEdges;

	while(ePrime != NullEdgeIdx)
	{
		wPrime = headOfEdgeId[ePrime];
		// eTupleUnWeighted eFlip(w, wPrime);
		// if(flippedEdges.find(ePrime) == flippedEdges.end())
		// {	
		flipDirectedEdge(ePrime, wPrime, w, EM); 	// flipDirctedEdge(eId, oldHeadNode, newHeadNode)
		w = wPrime;
		lastEId = ePrime;
		ePrime = getTightOutNbr(w);
		flippedEdges[ePrime] = 1;
		// }
		// else
		// {
		// 	break;
		// }
	}

	decrementDu(w, EM);
	// }
	return w;
}


VertexIdx DynamicGraph :: getMinDegreeVertexInE(edgeVector &currentEdge)
{
	VertexIdx minDegVertex = currentEdge[0];
	Count minDegree = nodeInDeg[currentEdge[0]];

	for(unsigned int i = 1; i < currentEdge.size(); i++)
	{
		if(nodeInDeg[currentEdge[i]] < minDegree)
		{
			minDegVertex = currentEdge[i];
			minDegree = nodeInDeg[currentEdge[i]];
		}
	}
	return minDegVertex;
}

int DynamicGraph :: addDirectedEdgeToInOutNbrs(edgeVector &e, EdgeIdx eId, VertexIdx newHeadNode)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	// VertexIdx headNode = newHeadNode;
	headOfEdgeId[eId] = newHeadNode;

	// add u to in-neighbors of v
	// InNbrs[v][u] += 1;
	InNbrs[newHeadNode].insert(eId);
	addEdgeToInNbrsForVisitNext(newHeadNode, eId);

	// add v to priority queue out-neighbors of u;
	// get the current value of d(v)
	Count headNodeVal = nodeInDeg[newHeadNode];

	// update e in the priority queue of all u's except v with this v's value...
	// edgeVector e = EM.edgeDupList[eId];
	for(unsigned int i = 0; i < e.size(); i++)
	{
		VertexIdx u = e[i];
		if(u != newHeadNode)
		{
			addToPriorityQueue(u, newHeadNode, headNodeVal, eId);
			InDegreeFromNodesView[u][newHeadNode] = headNodeVal;
		}
	}

	return 0;
}

std::pair<VertexIdx, EdgeIdx> DynamicGraph :: getTightInNbr(VertexIdx u, EdgeManager &EM)
{
	EdgeIdx neighborEdgeIdToReturn = NullEdgeIdx;
	VertexIdx neighborNodeIdToReturn = NullVertexIdx;

	std::list<EdgeIdx> :: iterator updateItInc;

	Count start = 0;
	Count maxNumNeighborsToCheck = (Count)(4 * nodeInDeg[u]) / eta;
	
	std::unordered_map<EdgeIdx, int> touchedNeighbors;
	// note that you shouldnt be updating an element multiple times within a same 
	// update... so keep track of the updated neighbors... 
	// every time start with a new/empty list of neighbors that would be updated...
	
	// this gives a pointer to a edgeId... 	
	updateItInc = nextPositionIteratorTightInNbr[u];

	while(start < maxNumNeighborsToCheck)
	{
		// access the next neighbor 
		EdgeIdx usNextNeighborEdgeId = *updateItInc;
		
		if(touchedNeighbors.find(usNextNeighborEdgeId) == touchedNeighbors.end())
		{
			// increase the counter...
			start++;
			// increment the pointer position...
			updateItInc++;

			// if we hit end of the list... round-robin to start of the list...
			if(updateItInc == listOfNeighbors[u].end())
	        {
	        	updateItInc = listOfNeighbors[u].begin();
	        }

			// add the updated neighbor to the list of updated neighbors...
			touchedNeighbors[usNextNeighborEdgeId] = 1;

			// check if this current neighbor satisfies the tight in-neighbor condition...
			// VertexIdx edgeHeadNode = headOfEdgeId[usNextNeighborEdgeId];
			edgeVector currentEdge = EM.edgeDupList[usNextNeighborEdgeId];
			VertexIdx minDegVertexInNbrE = getMinDegreeVertexInE(currentEdge);

			Count minNodeInDeg = nodeInDeg[minDegVertexInNbrE];
			Count newHeadNodeInDeg = nodeInDeg[u];
			
			if(minNodeInDeg <= newHeadNodeInDeg - eta/2)
			{
				neighborEdgeIdToReturn = usNextNeighborEdgeId;
				neighborNodeIdToReturn = minDegVertexInNbrE;
				break;
			}
		}
		else
		{
			// which says that we have already updated all the neigbors in the list...
			break;
		}
	}

	// update the position of the nextPosition Iterator... 
	nextPositionIteratorTightInNbr[u] = updateItInc;

	return std::make_pair(neighborNodeIdToReturn, neighborEdgeIdToReturn);
}

EdgeIdx DynamicGraph :: getTightOutNbr(VertexIdx u)
{
	EdgeIdx maxOutE = getMaxOutNbr(u);
	if(maxOutE != NullEdgeIdx)
	{
		VertexIdx t = headOfEdgeId[maxOutE];
		// degree of t in the view of u
		// if the max neighbor has the degree that is very high than that of u...
		float threshold = nodeInDeg[u] + eta/2;
		if((InDegreeFromNodesView[u][t] >= threshold))
		// if((nodeInDeg[t] >= nodeInDeg[u] + eta/2))
		{
			return maxOutE;
		}
	}
	return NullEdgeIdx;
}

EdgeIdx DynamicGraph :: getMaxOutNbr(VertexIdx u)
{
	// get the degToNode std::map of u
	if(outdegToNodeMap[u].size() >= 1)
	{
		std::map<Count, std::set<EdgeIdx>>::reverse_iterator rit = outdegToNodeMap[u].rbegin();
		Count maxVal = rit->first;
		std::set<EdgeIdx> maxValSet = rit->second;
		std::set<EdgeIdx>::iterator it = maxValSet.begin();
		EdgeIdx maxEle = *it;

		return maxEle;
	}

	return NullVertexIdx;
}


int DynamicGraph :: flipDirectedEdge(EdgeIdx eId, VertexIdx oldHeadNode, VertexIdx newHeadNode, EdgeManager &EM)
{
	edgeVector currentEdge = EM.edgeDupList[eId];
	removeDirectedEdgeFromInOutNbrs(currentEdge, eId, oldHeadNode);

	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	// eTupleUnWeighted flippedEdge (v,u);
	addDirectedEdgeToInOutNbrs(currentEdge, eId, newHeadNode);
	return 0;
}


int DynamicGraph :: incrementDu(VertexIdx headNode, EdgeManager &EM)
{
	updateLabels(headNode, 1);
	// Count oldVal = nodeInDeg[headNode];
	nodeInDeg[headNode] += 1;

	Count newDuVal = nodeInDeg[headNode];

	// update 4 din(headNode)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(headNode, newDuVal, 1, EM);

	return 0;
}

int DynamicGraph :: decrementDu(VertexIdx oldHeadNode, EdgeManager &EM)
{
	updateLabels(oldHeadNode, -1);
	// Count oldVal = nodeInDeg[oldHeadNode];
	nodeInDeg[oldHeadNode] -= 1;
	Count newDuVal = nodeInDeg[oldHeadNode];

	// update 4 din(u)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(oldHeadNode, newDuVal, -1, EM);
	
	return 0;
}


int DynamicGraph :: addEdgeToInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId)
{
	// This is to add v to the in-neighbors of u...
	// add to the std::list of in-neighbors and update the existence and 
	// address of this new neigbor in the map...
	listOfNeighbors[headNode].push_back(eId);
	std::list<EdgeIdx>::iterator esAddress = listOfNeighbors[headNode].end();
	--esAddress;									// going to the last element -- which is the new one inserted...
	mapToNeighborsList[headNode][eId] = esAddress;

	// if this is the first element that is getting added to the list...
	// let the next-pointers point to the first element/neighbor....
	if(listOfNeighbors[headNode].size() == 1)
	{
		nextPositionIteratorInc[headNode] = esAddress;
		nextPositionIteratorDec[headNode] = esAddress;
		nextPositionIteratorTightInNbr[headNode] = esAddress;
	}

	return 0;
}

int DynamicGraph :: addToPriorityQueue(VertexIdx u, VertexIdx headNode, Count headVal, EdgeIdx headEId)
{
	// we need to add/update headNode in the priority queue of u;
	if(nodeToOutdegMap[u].find(headEId) != nodeToOutdegMap[u].end())
	{
		Count oldVal = nodeToOutdegMap[u][headEId];
		nodeToOutdegMap[u][headEId] = headVal;

		outdegToNodeMap[u][oldVal].erase(headEId);
		if(outdegToNodeMap[u][oldVal].size() < 1)
		{
			outdegToNodeMap[u].erase(oldVal);
		}
		outdegToNodeMap[u][headVal].insert(headEId);
	}
	else
	{
		nodeToOutdegMap[u][headEId] = headVal;
		outdegToNodeMap[u][headVal].insert(headEId);
	}

	return 0;
}


int DynamicGraph :: removeDirectedEdgeFromInOutNbrs(edgeVector &e, EdgeIdx eId, VertexIdx oldHeadNode)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	// VertexIdx oldHeadNode = headOfEdgeId[eId];

	// remove/decrement u from in-neighbors of v
	// InNbrs[v][u] -= 1;

	InNbrs[oldHeadNode].erase(eId);
	removeEdgeFromInNbrsForVisitNext(oldHeadNode, eId);
	// edgeVector e = EM.edgeDupList[eId];

	for(unsigned int i = 0; i < e.size(); i++)
	{
		VertexIdx u = e[i];
		if(u != oldHeadNode)
		{
			removeFromPriorityQueue(u, oldHeadNode, eId);
		}
	}
	return 0;
}


int DynamicGraph :: updateLabels(VertexIdx u, Count changeVal)
{
	Count oldVal = nodeInDeg[u];
	Count newVal = nodeInDeg[u] + changeVal;
	// updating the Labels data-structure...
	Labels[oldVal].erase(u);
	// nodeInDeg[u] +=  changeVal;
	if(Labels[oldVal].size() == 0)
	{
		Labels.erase(oldVal);	
	}
	Labels[newVal].insert(u);
	ReverseLabels[u] = newVal;

	return 0;
}

int DynamicGraph :: updateNextNeighbors(VertexIdx headNode, Count newDuVal, int incOrDec, EdgeManager &EM)
{
	std::list<EdgeIdx> :: iterator updateItInc;
	std::unordered_map<EdgeIdx, int> touchedNeighbors;

	Count start = 0;
	Count maxNumNeighborsToUpdate = (Count) (4 * nodeInDeg[headNode]) / eta;
	// note that you shouldnt be updating an element multiple times within a same 
	// update... so keep track of the updated neighbors... 
	// every time start with a new/empty list of neighbors that would be updated...
	if(incOrDec == 1)
	{
		updateItInc = nextPositionIteratorInc[headNode];
	}
	else if (incOrDec == -1)
	{
		updateItInc = nextPositionIteratorDec[headNode];
	}

	if(updateItInc == listOfNeighbors[headNode].end())
    {
    	updateItInc = listOfNeighbors[headNode].begin();
    }

	while(start < maxNumNeighborsToUpdate)
	{
		// access the next neighbor 
		EdgeIdx usNextNeighbor = *updateItInc;
		
		if(touchedNeighbors.find(usNextNeighbor) == touchedNeighbors.end())
		{		
			// update new in-degree value of u to the neighbor
			// addToPriorityQueue(u, usNextNeighbor, newDuVal);
			edgeVector eNbr = EM.edgeDupList[usNextNeighbor];
			for(unsigned int i = 0; i < eNbr.size(); i++)
			{
				VertexIdx nbrNode = eNbr[i];
				if(nbrNode != headNode)
				{
					// addToPriorityQueue(nbrNode, headNode, newDuVal);
					InDegreeFromNodesView[nbrNode][headNode] = newDuVal;

					// remove eNbr from NbrNode 
					// add eNbr to NbrNode with new value... which is newDuVal;
					Count oldVal = nodeToOutdegMap[nbrNode][usNextNeighbor];
					outdegToNodeMap[nbrNode][oldVal].erase(usNextNeighbor);
					if(outdegToNodeMap[nbrNode][oldVal].size() < 1)
					{
						outdegToNodeMap[nbrNode].erase(oldVal);	
					}
					outdegToNodeMap[nbrNode][newDuVal].insert(usNextNeighbor);
					nodeToOutdegMap[nbrNode][usNextNeighbor] = newDuVal;
				}
			}
			// VertexIdx nbrHeadNode = headOfEdgeId[usNextNeighbor];

			// increase the counter...
			start++;
	
			// increment the pointer position...
			updateItInc++;
			// if we hit end of the list... round-robbin to start of the list...
			if(updateItInc == listOfNeighbors[headNode].end())
	        {
	        	updateItInc = listOfNeighbors[headNode].begin();
	        }

			// add the updated neighbor to the list of updated neighbors...
			touchedNeighbors[usNextNeighbor] = 1;
		}
		else
		{
			// which says that we have already updated all the neigbors in the list...
			break;
		}
	}

	// update the position of the nextPosition Iterator... 
	if(incOrDec == 1)
	{
		nextPositionIteratorInc[headNode] = updateItInc;
	}
	else
	{
		nextPositionIteratorDec[headNode] = updateItInc;
	}

	return 0;
}

int DynamicGraph :: removeEdgeFromInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId)
{
	// This is to remove v from the in-neighbors of u...

	// check if the current nextIterator is not at the same position as the 
	// node to be removed...
	// if so shift the position of the next-iterator to the next element
	if(mapToNeighborsList[headNode].find(eId) != mapToNeighborsList[headNode].end())
	{
		// if the address of the next position iterator matches with that of that of the 
		// element to be removed...
		if(nextPositionIteratorInc[headNode] == mapToNeighborsList[headNode][eId])
		{
			updateIncPointer(headNode);
		}

		if(nextPositionIteratorDec[headNode] == mapToNeighborsList[headNode][eId])
		{
			updateDecPointer(headNode);
		}

		if(nextPositionIteratorTightInNbr[headNode] == mapToNeighborsList[headNode][eId])
		{
			updateTightInNbrIterator(headNode);
		}

		// removing eId from the in-neighbors of u....
		// Remove eId using the iterator of eId, otherwise it is linear time... 
		// we want constant time...
		listOfNeighbors[headNode].erase(mapToNeighborsList[headNode][eId]);
		// remove element and its address from the map as well
		mapToNeighborsList[headNode].erase(eId);
	}

	return 0;
}

int DynamicGraph :: removeFromPriorityQueue(VertexIdx u, VertexIdx oldHeadNode, EdgeIdx oldHeadEId)
{
	Count oldVal = nodeToOutdegMap[u][oldHeadEId];
	outdegToNodeMap[u][oldVal].erase(oldHeadEId);
	if(outdegToNodeMap[u][oldVal].size() < 1)
	{
		outdegToNodeMap[u].erase(oldVal);
	}
	nodeToOutdegMap[u].erase(oldHeadEId);

	return 0;
}

int DynamicGraph :: updateIncPointer(VertexIdx u)
{
	// shift the itertor to the next position...
	std::list<VertexIdx>:: iterator updateItInc;

	updateItInc = nextPositionIteratorInc[u];
	updateItInc++;
	if(updateItInc == listOfNeighbors[u].end())
	{
		updateItInc = listOfNeighbors[u].begin();
	}
	nextPositionIteratorInc[u] = updateItInc;

	return 0;
}

int DynamicGraph :: updateDecPointer(VertexIdx u)
{
	// shift the itertor to the next position...
	std::list<VertexIdx>:: iterator updateItDec;

	updateItDec = nextPositionIteratorDec[u];
	updateItDec++;
	if(updateItDec == listOfNeighbors[u].end())
	{
		updateItDec = listOfNeighbors[u].begin();
	}
	nextPositionIteratorDec[u] = updateItDec;

	return 0;
}

int DynamicGraph :: updateTightInNbrIterator(VertexIdx u)
{
	// shift the itertor to the next position...
	std::list<VertexIdx>:: iterator updateItNext;

	updateItNext = nextPositionIteratorTightInNbr[u];
	updateItNext++;
	if(updateItNext == listOfNeighbors[u].end())
	{
		updateItNext = listOfNeighbors[u].begin();
	}
	nextPositionIteratorTightInNbr[u] = updateItNext;

	return 0;
}

double DynamicGraph :: getDensity()
{
	std::map<Count, std::set<VertexIdx>>::reverse_iterator labelsIt = Labels.rbegin();
	double currMaxDensityValue = (labelsIt->first)*(1-epsVal);
	// double currMaxDensityValue = (labelsIt->first);
	return currMaxDensityValue;
}

double DynamicGraph :: getRhoEst()
{
	return rhoEst;
}

Count DynamicGraph :: getMinLoadInE(EdgeIdx eId, EdgeManager &EM)
{
	edgeVector e = EM.edgeDupList[eId];
	VertexIdx minDegVertex = e[0];
	Count minDegree = nodeInDeg[e[0]];

	for(unsigned int i = 1; i < e.size(); i++)
	{
		if(nodeInDeg[e[i]] < minDegree)
		{
			minDegVertex = e[i];
			minDegree = nodeInDeg[e[i]];
		}
	}
	return minDegree;
}

int DynamicGraph :: addEdgeToPendingList(edgeVector &e, EdgeIdx eId)
{
	pendingListOfEdges[eId] = 1;

	for(unsigned int i = 0; i < e.size(); i++)
	{
		pendingForNode[e[i]].insert(eId);
		nodesWithPendingEdge[eId].insert(e[i]);
	}

	return 0;
}

int DynamicGraph :: checkEdgeExistenceInPendingList(EdgeIdx eId)
{
	if(pendingListOfEdges.find(eId) != pendingListOfEdges.end() && pendingListOfEdges[eId] == 1)  
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


// std::pair<EdgeIdx, edgeVector> Graph :: getPendingEdgeForLastVertex(VertexIdx lv)
EdgeIdx DynamicGraph :: getPendingEdgeForLastVertex(VertexIdx lv)
{
	// edgeVector e = {};
	EdgeIdx eId = NullEdgeIdx;

	// std::pair<EdgeIdx, edgeVector> eIdEdgePair = std::make_pair(eId, e);

	std::unordered_map<VertexIdx, std::set<VertexIdx>> :: iterator pit = pendingForNode.find(lv);

	if(pit != pendingForNode.end())
	{
		std::set<EdgeIdx> edgeSet = pit->second;
		std::set<EdgeIdx> :: iterator setIt = edgeSet.begin();
		eId = *setIt;
		// e = pendingListOfEdges[eId];
		// eIdEdgePair = std::make_pair(eId, e);
	}

	// return eIdEdgePair;
	return eId;
}


int DynamicGraph :: removeEdgeFromPendingList(edgeVector &e, EdgeIdx eId)
{
	for(unsigned int i = 0; i < e.size(); i++)
	{
		pendingForNode[e[i]].erase(eId);
		nodesWithPendingEdge[eId].erase(e[i]);
	}

	pendingListOfEdges.erase(eId);

	return 0;
}