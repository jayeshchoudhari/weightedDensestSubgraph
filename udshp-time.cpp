#include "Graph.h"
#include <bits/stdc++.h>
#include <chrono>

using namespace rwNameSpace;

class Graph 
{
    private:
        static Count nVertices; 								//number of vertices in the graph
        static Count nEdges;     								//number of edges in this list
        static double alpha;										// multiplicity
        static double epsVal;									// error factor

        static std::multimap<std::vector<VertexIdx>, EdgeIdx> edgeMap;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        static std::vector <std::vector<VertexIdx>> edgeList;

        // above could be static members shared among all the instances of the Graph...
        // each instance for a guess of density...
	
        // std::multimap<std::vector<VertexIdx>, EdgeIdx> pendingListOfEdges;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        std::unordered_map<EdgeIdx, std::vector<VertexIdx>> pendingListOfEdges;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        std::unordered_map<VertexIdx, std::set<EdgeIdx>> pendingForNode;
        std::unordered_map<EdgeIdx, std::set<VertexIdx>> nodesWithPendingEdge;

 		double eta;
 		double rhoEst;

        std::unordered_map<EdgeIdx, VertexIdx> headOfEdgeId;

        std::map <VertexIdx, Count> nodeInDeg;

        std::vector <std::set<EdgeIdx>> InNbrs;							//List of InNbrs

        // ****** FOR MAINTAINING THE VISITNEXT DATA-STRUCTURE ********
        // list of neighbors for each node...
        // vector<std::list<VertexIdx>> listOfNeighbors;
        std::vector<std::list<EdgeIdx>> listOfNeighbors;
        // map to the elements in the neighbors list for each node...
        std::vector<std::map<EdgeIdx, std::list<VertexIdx>::iterator>> mapToNeighborsList;
        // std::vector of iterators -- one for each node... that keeps track of which element to access next in visitNext
        std::vector<std::list<EdgeIdx>::iterator> nextPositionIteratorInc;
        std::vector<std::list<EdgeIdx>::iterator> nextPositionIteratorDec;
        std::vector<std::list<EdgeIdx>::iterator> nextPositionIteratorTightInNbr;
        // *************************************************************

        // ********* TO MAINTAIN A UPDATING PRIORITY QUEUE OF OUTNEIGHBORS *************** 
        // this is to be done for each node.. and thus a vector here...
        // Each element of a vector is a map --  where the map is with a key Count, and value is a set of 
        // of edgeIds, which have the headVertices with value (indegree) as Count...
        // and nodeToOutdegMap is a reverseMap of the same...
        std::vector <std::map<Count, std::set<EdgeIdx>> > outdegToNodeMap;
        std::vector <std::map<EdgeIdx, Count>> nodeToOutdegMap;
        std::vector <std::unordered_map<VertexIdx, Count>> InDegreeFromNodesView;
        // ********************************************************************************

        // ********** TO MAINTAIN LABELS/NODES WITH HIGH INDEGREE *************** 
        std::map<Count, std::set<VertexIdx>> Labels;				// this is the data-structure that keeps track of top elements...
        std::map<VertexIdx, Count> ReverseLabels;			// this is the data-structure that keeps track of top elements...
        // **********************************************************************

    public:

    	static int initializeStaticMembers(Count nv, double decEps);

    	int initializeVariables(int i, Count nv);

    	int showInstanceVariables();

    	double getRhoEst();

    	int addEdgeToPendingList(edgeVector e, EdgeIdx eId);
    	int checkEdgeExistenceInPendingList(edgeVector e, EdgeIdx eId);
		int removeEdgeFromPendingList(edgeVector e, EdgeIdx eId);
		std::pair<EdgeIdx, edgeVector> getPendingEdgeForLastVertex(VertexIdx lv);

		// Graph(int nv, int ne);		//gets number of vertices and edges... edges might change -- but this is just about the file... 

		// Graph(Count nv, double decEta, double decEps);		
		//gets number of vertices and edges... edges might change -- but this is just about the file... 
    	// Graph(std::vector <VertexIdx> &v, std::vector<edgeVector> &edList);		//takes in list/std::vector of nodes and std::vector of pairs which are edges (undirected)

		Count getNumVertices();
		Count getNumEdges();
        // //Make a copy
        // Graph copy() const;

        // //For debugging.  Not for serialization.
        // void print(FILE *f = stdout) const;

		int addDirectedEdgeToInOutNbrs(EdgeIdx eId, VertexIdx v);
		int removeDirectedEdgeFromInOutNbrs(EdgeIdx eId, VertexIdx headNode);
		int flipDirectedEdge(EdgeIdx eId, VertexIdx oldHeadNode, VertexIdx newHeadNode);
		
		int addToPriorityQueue(VertexIdx, VertexIdx, Count, EdgeIdx);
		int removeFromPriorityQueue(VertexIdx, VertexIdx, EdgeIdx);
		
		int updateNextNeighbors(VertexIdx u, Count newDuVal, int incOrDec);

		int incrementDu(VertexIdx);
		int decrementDu(VertexIdx);

		std::pair<VertexIdx, EdgeIdx> getTightInNbr(VertexIdx);
		VertexIdx getMinDegreeVertexInE(EdgeIdx eId);
		Count getMinLoadInE(EdgeIdx eId);
		EdgeIdx getTightOutNbr(VertexIdx);
		EdgeIdx getMaxOutNbr(VertexIdx u);

		int updateLabels(VertexIdx u, Count changeVal);

		int insertEdge(edgeVector e, EdgeIdx eId);
		int deleteEdge(edgeVector e, EdgeIdx eId);
		VertexIdx deleteEdgeReturnLastVertex(edgeVector e, EdgeIdx eId);

		int updateIncPointer(VertexIdx u);
		int updateDecPointer(VertexIdx u);
		int updateTightInNbrIterator(VertexIdx u);
		int addEdgeToInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId);
		int removeEdgeFromInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId);

		int setEdgeId(edgeVector e, EdgeIdx eId);
		EdgeIdx checkEdgeExistence(edgeVector e);
		EdgeIdx getEdgeId(edgeVector e);
		int removeEdgeFromMap(edgeVector e);
		int addEdgeToEdgeList(edgeVector e);

		Count getLabel(VertexIdx);
		Count getMaxLabel();
		double getDensity();

		std::pair<std::set<VertexIdx>, revItMapCountSetVertices> returnDensitySatisfiedNodes(revItMapCountSetVertices startIt, Count D);
		std::set<VertexIdx> getDensestSubgraph(double rgamma);
		std::set<VertexIdx> querySubgraph(double D_hat);
		Count getMaxIndegree();
		int checkEdgeAssignment();

		double getDensityOfInducedSubgraph(std::set<VertexIdx>);
		std::pair<double, unsigned int> getMaxPartitionDensity();

		unsigned int getPendingCount();
		int addListOfPendingEdges();

		int showPQs();
};


// initializing static (class) members...
Count Graph::nVertices{0};
Count Graph::nEdges{0};
double Graph :: alpha{0};
double Graph :: epsVal{0};
std::multimap<std::vector<VertexIdx>, EdgeIdx> Graph::edgeMap;
std::vector <std::vector<VertexIdx>> Graph::edgeList;


// Constructor to initialize the graph...
// Graph :: Graph(int nv, int ne)

int Graph :: initializeStaticMembers(Count nv, double decEps)
{
    nVertices = nv;
	epsVal = decEps;
    alpha = round(1.0/(pow(epsVal,2)));

    return 0;
}

int Graph :: initializeVariables(int i, Count nv)
{

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

	outdegToNodeMap.resize(nv);
	nodeToOutdegMap.resize(nv);
	InDegreeFromNodesView.resize(nv);

	InNbrs.resize(nv);

	// to maintain visitNext structure... 
	listOfNeighbors.resize(nv);
	mapToNeighborsList.resize(nv);
	nextPositionIteratorInc.resize(nv);
	nextPositionIteratorDec.resize(nv);
	nextPositionIteratorTightInNbr.resize(nv);

	return 0;
}


int Graph :: showPQs()
{
	VertexIdx ni = 0;
	std::cout << "Printing outdegToNodeMap -- std::vector <std::map<Count, std::set<EdgeIdx>>>\n"; 
	std::vector <std::map<Count, std::set<EdgeIdx>>> :: iterator it1;
	for(it1 = outdegToNodeMap.begin(); it1 != outdegToNodeMap.end(); it1++)
	{
		std::cout << ni << "\n";
		std::map<Count, std::set<EdgeIdx>> temp = *it1;
		std::map<Count, std::set<EdgeIdx>>::iterator inIt1;
		for(inIt1 = temp.begin(); inIt1 != temp.end(); inIt1++)
		{
			std::cout << "---" << inIt1->first << " : "; 
			std::set<EdgeIdx> tset = inIt1->second;
			std::set<EdgeIdx>::iterator tsetIt;
			for(tsetIt = tset.begin(); tsetIt != tset.end(); tsetIt++)
			{
				std::cout << *tsetIt << ", ";
			}
			std::cout << std::endl;
		}
		ni++;
	}

	std::cout << "Printing nodeToOutdegMap -- vector <std::map<EdgeIdx, Count>>\n";
	ni = 0;
    std::vector <std::map<EdgeIdx, Count>>::iterator it2;
    for(it2 = nodeToOutdegMap.begin(); it2 != nodeToOutdegMap.end(); it2++)
    {
    	std::cout << ni << "\n";
    	std::map<EdgeIdx, Count> temp = *it2;
    	std::map<EdgeIdx, Count>::iterator inIt1;
    	for(inIt1 = temp.begin(); inIt1 != temp.end(); inIt1++)
    	{
    		std::cout << " --- " << inIt1->first << " " << inIt1->second << "\n"; 
    	}

    	ni++;
    }


	return 0;
}

int Graph :: addEdgeToEdgeList(edgeVector e)
{
	edgeList.push_back(e);
	return 0;
} 

int Graph :: addEdgeToPendingList(edgeVector e, EdgeIdx eId)
{
	pendingListOfEdges[eId] = e;

	for(unsigned int i = 0; i < e.size(); i++)
	{
		pendingForNode[e[i]].insert(eId);
		nodesWithPendingEdge[eId].insert(e[i]);
	}

	return 0;
}

int Graph :: checkEdgeExistenceInPendingList(edgeVector e, EdgeIdx eId)
{
	if(pendingListOfEdges.find(eId) != pendingListOfEdges.end())
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


std::pair<EdgeIdx, edgeVector> Graph :: getPendingEdgeForLastVertex(VertexIdx lv)
{
	edgeVector e = {};
	EdgeIdx eId = NullEdgeIdx;

	std::pair<EdgeIdx, edgeVector> eIdEdgePair = std::make_pair(eId, e);

	std::unordered_map<VertexIdx, std::set<VertexIdx>> :: iterator pit = pendingForNode.find(lv);

	if(pit != pendingForNode.end())
	{
		std::set<EdgeIdx> edgeSet = pit->second;
		std::set<EdgeIdx> :: iterator setIt = edgeSet.begin();
		eId = *setIt;
		e = pendingListOfEdges[eId];
		eIdEdgePair = std::make_pair(eId, e);
	}

	return eIdEdgePair;
}


int Graph :: removeEdgeFromPendingList(edgeVector e, EdgeIdx eId)
{
	for(unsigned int i = 0; i < e.size(); i++)
	{
		pendingForNode[e[i]].erase(eId);
		nodesWithPendingEdge[eId].erase(e[i]);
	}

	pendingListOfEdges.erase(eId);
	
	// std::multimap<std::vector<VertexIdx>, EdgeIdx>:: iterator pit = pendingListOfEdges.find(e);
	// if(pit != pendingListOfEdges.end())
	// {
	// 	pendingListOfEdges.erase(pit);
	// }
	return 0;
}




Count Graph :: getNumVertices()
{
	return nVertices;
}

int Graph :: showInstanceVariables()
{
	std::cout << "num nodes = " << nVertices << "\n";
    std::cout << "eta = " << eta << "\n";
    std::cout << "epsVal = " << epsVal << "\n";
    std::cout << "rhoEst = " << rhoEst << "\n";
    std::cout << "Alpha = " << alpha << "\n";

	return 0;
}


EdgeIdx Graph :: checkEdgeExistence(edgeVector e)
{
	std::multimap<std::vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
	if(edgeMapIt == edgeMap.end())
	{
		return NullEdgeIdx;
	}
	else
	{
		EdgeIdx eId = edgeMapIt->second;
		return eId;
	}
}


int Graph :: setEdgeId(edgeVector e, EdgeIdx eId)
{
	// edgeMap[e] = eId;
	edgeMap.insert(std::pair<std::vector<VertexIdx>, EdgeIdx>(e, eId));
	// reverseEdgeMap[eId] = e;
	return 0;
}

EdgeIdx Graph :: getEdgeId(edgeVector e)
{
	std::multimap<std::vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
	EdgeIdx eId = edgeMapIt->second;
	return eId;
}

int Graph :: removeEdgeFromMap(edgeVector e)
{
	// edgeMap.erase(e);
	std::multimap<std::vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
	edgeMap.erase(edgeMapIt);
	return 0;
}


int Graph :: updateNextNeighbors(VertexIdx headNode, Count newDuVal, int incOrDec)
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
			edgeVector eNbr = edgeList[usNextNeighbor];
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

int Graph :: updateLabels(VertexIdx u, Count changeVal)
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

int Graph :: incrementDu(VertexIdx headNode)
{
	updateLabels(headNode, 1);
	// Count oldVal = nodeInDeg[headNode];
	nodeInDeg[headNode] += 1;

	Count newDuVal = nodeInDeg[headNode];

	// update 4 din(headNode)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(headNode, newDuVal, 1);

	return 0;
}

int Graph :: decrementDu(VertexIdx oldHeadNode)
{
	updateLabels(oldHeadNode, -1);
	// Count oldVal = nodeInDeg[oldHeadNode];
	nodeInDeg[oldHeadNode] -= 1;
	Count newDuVal = nodeInDeg[oldHeadNode];

	// update 4 din(u)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(oldHeadNode, newDuVal, -1);
	
	return 0;
}

int Graph :: addToPriorityQueue(VertexIdx u, VertexIdx headNode, Count headVal, EdgeIdx headEId)
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

// int Graph :: addToInNbrs(VertexIdx u, VertexIdx v)
int Graph :: addEdgeToInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId)
{
	// This is to add v to the in-neighbors of u...
	// add to the std::list of in-neighbors and update the existence and 
	// address of this new neigbor in the map...
	listOfNeighbors[headNode].push_back(eId);				// adding v to the in-neighbors of u...
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

int Graph :: updateIncPointer(VertexIdx u)
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

int Graph :: updateDecPointer(VertexIdx u)
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


int Graph :: updateTightInNbrIterator(VertexIdx u)
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

// int Graph :: removeFromInNbrs(VertexIdx u, VertexIdx v)
int Graph :: removeEdgeFromInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId)
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


int Graph :: addDirectedEdgeToInOutNbrs(EdgeIdx eId, VertexIdx newHeadNode)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	VertexIdx headNode = newHeadNode;
	headOfEdgeId[eId] = headNode;

	// add u to in-neighbors of v
	// InNbrs[v][u] += 1;
	InNbrs[headNode].insert(eId);
	addEdgeToInNbrsForVisitNext(headNode, eId);

	/*
	// Note that this is to be done only for the first time when u becomes in-neighbor of v..
	// remember this is a multigraph.... 
	if(InNbrs[v][u] == 1)
	{
		addToInNbrs(v, u);
	}
	*/

	// add v to priority queue out-neighbors of u;
	// get the current value of d(v)
	// Count headNodeVal = nodeInDeg[headNode] + 1;
	Count headNodeVal = nodeInDeg[headNode];

	// update e in the priority queue of all u's except v with this v's value...
	edgeVector e = edgeList[eId];
	for(unsigned int i = 0; i < e.size(); i++)
	{
		VertexIdx u = e[i];
		if(u != headNode)
		{
			addToPriorityQueue(u, headNode, headNodeVal, eId);
			InDegreeFromNodesView[u][headNode] = headNodeVal;
		}
	}

	return 0;
}

int Graph :: removeFromPriorityQueue(VertexIdx u, VertexIdx oldHeadNode, EdgeIdx oldHeadEId)
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


int Graph :: removeDirectedEdgeFromInOutNbrs(EdgeIdx eId, VertexIdx oldHeadNode)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	// VertexIdx oldHeadNode = headOfEdgeId[eId];

	// remove/decrement u from in-neighbors of v
	// InNbrs[v][u] -= 1;

	InNbrs[oldHeadNode].erase(eId);
	removeEdgeFromInNbrsForVisitNext(oldHeadNode, eId);
	edgeVector e = edgeList[eId];

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

int Graph :: flipDirectedEdge(EdgeIdx eId, VertexIdx oldHeadNode, VertexIdx newHeadNode)
{
	removeDirectedEdgeFromInOutNbrs(eId, oldHeadNode);

	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	// eTupleUnWeighted flippedEdge (v,u);
	addDirectedEdgeToInOutNbrs(eId, newHeadNode);
	return 0;
}

VertexIdx Graph :: getMinDegreeVertexInE(EdgeIdx eId)
{
	edgeVector e = edgeList[eId];
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
	return minDegVertex;
}

rwNameSpace::Count Graph :: getMinLoadInE(EdgeIdx eId)
{
	edgeVector e = edgeList[eId];
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


std::pair<VertexIdx, EdgeIdx> Graph :: getTightInNbr(VertexIdx u)
{
	// VertexIdx neighborToReturn = NullVertexIdx;
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

			// if we hit end of the list... round-robbin to start of the list...
			if(updateItInc == listOfNeighbors[u].end())
	        {
	        	updateItInc = listOfNeighbors[u].begin();
	        }

			// add the updated neighbor to the list of updated neighbors...
			touchedNeighbors[usNextNeighborEdgeId] = 1;

			// check if this current neighbor satisfies the tight in-neighbor condition...
			// VertexIdx edgeHeadNode = headOfEdgeId[usNextNeighborEdgeId];
			VertexIdx minDegVertexInNbrE = getMinDegreeVertexInE(usNextNeighborEdgeId);

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

EdgeIdx Graph :: getTightOutNbr(VertexIdx u)
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



EdgeIdx Graph :: getMaxOutNbr(VertexIdx u)
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



int Graph :: insertEdge(edgeVector e, EdgeIdx eId)
{
	// e = u,v -- at this moment an edge isnt directed...
	// hereafter it will be... from here onwards we would orient the edges...
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	VertexIdx w, wPrime;
	EdgeIdx ePrime, lastEId;

	w = getMinDegreeVertexInE(eId);
	lastEId = eId;

	// std::cout << nodeInDeg[w] << " headNode indegree during addition\n";
	addDirectedEdgeToInOutNbrs(eId, w);

	// std::cout << InNbrs[w].size() << " headNode innbrs size during addition\n";

	// headOfEdgeId[eId] = w;

	// check if this results into making some neighboring edge of w tight...
	std::pair<VertexIdx, EdgeIdx> minNbrNodeEdgePair = getTightInNbr(w);
	ePrime = minNbrNodeEdgePair.second;		// tight edge... -- head of this edge give tight In-neighbor node...
	// i.e. an edge with a headNode whose degree is very less as compared to that of w....
	// while(wPrime != -1)
	while(ePrime != -1)
	{
		wPrime = minNbrNodeEdgePair.first;
		// wPrime's in-deg is very less...
		// so flip the edge to wPrime..?

		// eTupleUnWeighted eFlip(wPrime, w);
		flipDirectedEdge(ePrime, w, wPrime);		// flipDirectedEdge(eId, oldHeadNode, newHeadNode)
		// so now wPrime is settled -- i.e. we have increased its indegree 
		// now we need to check if any in-neighbor of wPrime, violates the condition or has very less indegree...
		w = wPrime;
		lastEId = ePrime;
		minNbrNodeEdgePair = getTightInNbr(w);
		ePrime = minNbrNodeEdgePair.second;
	}

	incrementDu(w);

	return 0;
}


int Graph :: deleteEdge(edgeVector e, EdgeIdx eId)
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
		removeDirectedEdgeFromInOutNbrs(eId, headNode);

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
			flipDirectedEdge(ePrime, wPrime, w); 	// flipDirctedEdge(eId, oldHeadNode, newHeadNode)
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

		decrementDu(w);
	// }

	return 0;
}


VertexIdx Graph :: deleteEdgeReturnLastVertex(rwNameSpace::edgeVector e, EdgeIdx eId)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	VertexIdx w, wPrime;
	EdgeIdx ePrime, lastEId;

	VertexIdx headNode;

	if(headOfEdgeId.find(eId) != headOfEdgeId.end())
	{
		headNode = headOfEdgeId[eId];
	}
	else
	{
		headNode = NullVertexIdx;
	}

	// std::cout << nodeInDeg[headNode] << " headNode indegree before deletion\n";
	// std::cout << InNbrs[headNode].size() << " headNode innbrs size before deletion\n";
	
	if(headNode != NullVertexIdx)
	{
		// remove e from the InNbrs of headNode...
		removeDirectedEdgeFromInOutNbrs(eId, headNode);

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
			flipDirectedEdge(ePrime, wPrime, w); 	// flipDirctedEdge(eId, oldHeadNode, newHeadNode)
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

		decrementDu(w);
	}

	return w;
}


double Graph :: getDensityOfInducedSubgraph(std::set<VertexIdx> denseSubgraphNodes)
{

	// std::cout << "Checking density of said nodes... -- ";
	Count i = 0;
	std::multimap<std::vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt;

	std::vector<edgeVector> selectedEdges;

	Count edgeCount = 0;

	std::vector<VertexIdx> evector;

	for(edgeMapIt = edgeMap.begin(); edgeMapIt != edgeMap.end(); ++edgeMapIt)
	{
		i++;
		evector = edgeMapIt->first;
		// EdgeIdx eId = edgeMapIt->second;
		int flag = 1;
		for(unsigned int i = 0; i < evector.size(); i++)
		{
			VertexIdx node = evector[i];
			if(denseSubgraphNodes.find(node) == denseSubgraphNodes.end())
			{
				flag = 0;
				break;
			}
		}

		if(flag == 1)
		{
			edgeCount += 1;
			selectedEdges.push_back(evector);
		}
	}

	double estDensity = (edgeCount*1.0)/denseSubgraphNodes.size();

	std::cout << i << "---" << (edgeCount*1.0)/denseSubgraphNodes.size() << " " << denseSubgraphNodes.size() << std::endl;

	return estDensity;
}


std::pair<double, unsigned int> Graph :: getMaxPartitionDensity()
{
	double maxDensity = -1;
	unsigned int maxSubgraphSize = 0;

	std::set<VertexIdx> denseSubgraphNodes;

	revItMapCountSetVertices rit;

	for(rit = Labels.rbegin(); rit != Labels.rend(); ++rit)
	{
		std::set<VertexIdx> newElementsToSet = rit->second;

		std::set<VertexIdx> :: iterator bIt;
		for(bIt = newElementsToSet.begin(); bIt != newElementsToSet.end(); ++bIt)
		{
			denseSubgraphNodes.insert(*bIt);
		}
		
		double currentDensity = getDensityOfInducedSubgraph(denseSubgraphNodes);
		
		if(currentDensity > maxDensity)
		{
			maxDensity = currentDensity;
			maxSubgraphSize = denseSubgraphNodes.size();
		}
	}

	std::cout << "pending list of edges size ---- " << pendingListOfEdges.size() << std::endl;

	return std::make_pair(maxDensity, maxSubgraphSize);
}


unsigned int Graph :: getPendingCount()
{
	return pendingListOfEdges.size();	
}

int Graph :: addListOfPendingEdges()
{

	std::cout << "Adding pending list of edges to active copy...\n";
	std::unordered_map<EdgeIdx, std::vector<VertexIdx>> :: iterator pendingIt;

	for(pendingIt = pendingListOfEdges.begin(); pendingIt != pendingListOfEdges.end(); ++pendingIt)
	{
		EdgeIdx pEId = pendingIt->first;
		edgeVector pE = pendingIt->second;

		insertEdge(pE, pEId);
	}

	pendingListOfEdges.clear();
	pendingForNode.clear();
	nodesWithPendingEdge.clear();

	return 0;
}

std::pair<std::set<VertexIdx>, revItMapCountSetVertices> Graph :: returnDensitySatisfiedNodes(revItMapCountSetVertices startIt, Count D)
{
	std::set<VertexIdx> B;
    std::map<Count, std::set<VertexIdx>>::reverse_iterator preservedRit;
    std::map<Count, std::set<VertexIdx>>::reverse_iterator rit;

	for(rit = startIt; rit != Labels.rend(); ++rit)
	{
		preservedRit = rit;
		Count densityVal = rit->first;
		if(densityVal >= D)
		{
			std::set<VertexIdx> dvertices = rit->second;
			std::set<VertexIdx> :: iterator dvIt;
			for(dvIt = dvertices.begin(); dvIt != dvertices.end(); ++dvIt)
			{
				B.insert(*dvIt);
			}
		}
		else
		{
			break;
		}
	}

	if(rit == Labels.rend())
	{
		preservedRit = rit;
	}

	std::pair rP = std::make_pair(B, preservedRit);
    return rP;
}



double Graph :: getDensity()
{
	std::map<Count, std::set<VertexIdx>>::reverse_iterator labelsIt = Labels.rbegin();
	double currMaxDensityValue = (labelsIt->first)*(1-epsVal);
	// double currMaxDensityValue = (labelsIt->first);
	return currMaxDensityValue;
}

std::set<VertexIdx> Graph :: getDensestSubgraph(double rgamma)
{
	Count decrementVal = eta;

	double D = getMaxLabel();
	std::set<VertexIdx> A, B;
	unsigned int ASize, BSize;
	double sizeRatio;

	revItMapCountSetVertices rit;

	rit = Labels.rbegin();
	std::pair<std::set<VertexIdx>, revItMapCountSetVertices> AElementsPair = returnDensitySatisfiedNodes(rit, D);
	A = AElementsPair.first;

	if(AElementsPair.second != Labels.rend())
	{
		rit = AElementsPair.second;
	}
	else
	{
		std::cout << "Returning just the elements in A -- which is a set of all elements...\n";
		return A;
	}

	// D = D-eta;
	D = D-decrementVal;
	B = A;

	std::pair<std::set<VertexIdx>, revItMapCountSetVertices> newElementsToBPair = returnDensitySatisfiedNodes(rit, D);
	std::set<VertexIdx> newElementsToB = newElementsToBPair.first;

	rit = newElementsToBPair.second;

	std::set<VertexIdx> :: iterator bIt;	

	for(bIt = newElementsToB.begin(); bIt != newElementsToB.end(); ++bIt)
	{
		B.insert(*bIt);
	}
	
	// for(unsigned int i = 0; i < newElementsToB.size(); i++)
	// {
	// 	B.insert(newElementsToB[i]);
	// }

	BSize = B.size();
	ASize = A.size();
	sizeRatio = (BSize * 1.0)/ASize;

	while(sizeRatio > 1 + rgamma)
	{
		if(rit != Labels.rend())
		{
			A = B;
			// D = D - eta;
			D = D - decrementVal;
			newElementsToBPair = returnDensitySatisfiedNodes(rit, D);
			newElementsToB = newElementsToBPair.first;
			rit = newElementsToBPair.second;

			// for(unsigned int i = 0; i < newElementsToB.size(); i++)
			// {
			// 	B.insert(newElementsToB[i]);
			// }

			for(bIt = newElementsToB.begin(); bIt != newElementsToB.end(); ++bIt)
			{
				B.insert(*bIt);
			}

			BSize = B.size();
			ASize = A.size();
			sizeRatio = (BSize * 1.0)/ASize;
		}
		else
		{
			std::cout << "Cannot accumulate B anymore...\n";
			return B;
		}
	}

	return B;
}


std::set<VertexIdx> Graph :: querySubgraph(double D_hat)
{
	double rgamma = sqrt(2 * eta * log(nVertices) / rhoEst);
	std::cout << "rGamma = " << rgamma << std::endl;
	return getDensestSubgraph(rgamma);
}


Count Graph :: getMaxLabel()
{
	std::map<Count, std::set<VertexIdx>> :: reverse_iterator rit = Labels.rbegin();
	Count maxVal = rit->first;
	/*
	std::set<VertexIdx> maxValSet = rit->second;
	std::set<VertexIdx>::iterator it = maxValSet.begin();
	VertexIdx maxEle = *it;
	*/
	return maxVal;
}


double Graph :: getRhoEst()
{
	return rhoEst;
}

Count Graph :: getLabel(VertexIdx u)
{
	return ReverseLabels[u];
}

// utilities...
int sampleFromBinomial(int wt, double p)
{
	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// default_random_engine generator (seed);
  	std::random_device rd;
    std::mt19937 gen(rd());
    // perform 4 trials, each succeeds 1 in 2 times
    // std::binomial_distribution<> d(4, 0.5);

  	std::binomial_distribution<int> distribution (wt, p);
  	return distribution(gen);
}


int printEdgeVector(edgeVector e)
{
	for(unsigned int i = 0; i < e.size(); i++)
	{
		std::cout << e[i] << " ";
	}
	std::cout << "\n";

	return 0;
}


int main(int argc, char** argv)
{
	std::string graphFileName = argv[1];
	std::string outFileName = argv[2];
	double epsUD = std::stod(argv[3]);

	double localAlpha = round(1.0/pow(epsUD,2));
	double rho, rhoEstActive, rhoEst;

	std::ifstream graphFile;
	std::ofstream outFile;

	std::string line;
	std::stringstream ss;
	
	Count n, m, maxKEdge, totalProcessedEdges = 0; 
	VertexIdx e_src, e_dest, eEle;
	Weight e_weight, Wmax;
	char insDel;

	double scalingProbParam_c, sparsifyProbParam_c, scalingEps, sparsifyEps;
	double scalingEpsSquared = scalingEps * scalingEps;
	double sparsifyEpsSquared = sparsifyEps * sparsifyEps;

	double scalingProb, minWeightedDensity;

	std::chrono::steady_clock::time_point beginClock;
	std::chrono::steady_clock::time_point endClock;

	long long int totalTime = 0;
	long long int timeDiff = 0;
	long long int localInsertTime = 0;
	long long int localDeleteTime = 0;
	long long int perReportTime = 0;

	int active = 0;

	// std::cout << "Did I come here...\n";

	graphFile.open(graphFileName, std::ifstream::in);
	if(graphFile.is_open())
	{
		// if (getline(graphFile, line) != -1)
		if (getline(graphFile, line))
		{
			outFile.open(outFileName, std::ofstream::out);

			outFile << "YearLabel" << " " << "MaxInDeg" << " " << "(1-eps)MaxInDeg" << " " << "EstimatedDensity" << " " << "DenseSubgraphSize" << " " << "MaxPartitionDensity" << " " << "MaxPartitionSize" << " " << "TotalProcessedEdges" << " " << "EdgeAdditions" << std::endl;

			ss.clear();
			ss.str("");
			ss << line;
			// ss >> n >> m >> Wmax;
			ss >> n >> maxKEdge;


			int logn = (int)log2(n);
			int numInstances = logn + 1;
			int maxInstanceId = numInstances - 1;

			// for now, let us say that only the exact number of nodes is only known... 
			// initialize the datastructure required for each node...
			// d(u) 
			// InNbrs -- List of in-nbrs for each node
			// OutNbrs -- Max priority queue for each node u... indexed by d(v)

			Graph Ti[numInstances];
			Ti[0].initializeStaticMembers(n, epsUD);

			for(int i = 0; i <= maxInstanceId; ++i)
			{
				Ti[i].initializeVariables(i, n);
				// Ti[i].showInstanceVariables();
			}

			// Graph G(n, etaVal, epsUD); 	// This initializes du, InNbrs, OutNbrs(list of priority queues...)
			std::cout << "Done with initialization...\n";
		
			minWeightedDensity = (Wmax*1.0)/2;

			// scalingProb = (c log n) / (eps^2 * rho_min)  
			scalingProb = (scalingProbParam_c * log2(n)) / (scalingEpsSquared * minWeightedDensity);

			EdgeIdx edgeId = 0;
			Count edgeAdditions = 0, edgeDeletions = 0;
			edgeVector currentEdge;

			while(getline(graphFile, line))
			{
				ss.clear();
				ss.str("");
				ss << line;

				ss >> insDel;

				edgeVector tempVec;

				if(insDel == '+')
				{
					while(ss >> eEle)
					{
						tempVec.push_back(eEle);
					}

					if(tempVec.size() >= 2)
					{
						//pop_back year val... the last index....
						Count yearVal = tempVec[tempVec.size()-1]; 	// last element is the year/timestamp value...
						tempVec.pop_back();

						currentEdge = tempVec;

						if(currentEdge.size() >= 2)
						{
							beginClock = std::chrono::steady_clock::now();
							for(int dup = 0; dup < localAlpha; ++dup)
							{
								edgeAdditions += 1;
								totalProcessedEdges += 1;

								Ti[0].setEdgeId(currentEdge, edgeId);
								Ti[0].addEdgeToEdgeList(currentEdge);
								
								// if(totalProcessedEdges > 20700000)
								// {
								// 	std::cout << "Add the edge -- " << edgeId << " -- "; 
								// 	printEdgeVector(currentEdge);
								// }

								for(int copy = maxInstanceId; copy >= active + 1; --copy)
								{
									Ti[copy].insertEdge(currentEdge, edgeId);
								}

								if(active + 1 <= maxInstanceId)
								{
									rho = Ti[active + 1].getDensity();
									rhoEstActive = Ti[active].getRhoEst();

									if(rho >= 2*rhoEstActive)
									{
										active = active + 1;
									}
									else
									{
										Ti[active].insertEdge(currentEdge, edgeId);
									}
								}
								else
								{
									Ti[active].insertEdge(currentEdge, edgeId);
								}

								for(int copy = active - 1; copy >= 1; --copy)
								{
									// Ti[copy].insertEdge(currentEdge, edgeId);
									
									rhoEst = Ti[copy].getRhoEst();
									Count minLoad = Ti[copy].getMinLoadInE(edgeId);

									if(minLoad >= 2 * rhoEst)
									{
										Ti[copy].addEdgeToPendingList(currentEdge, edgeId);
									}
									else
									{
										Ti[copy].insertEdge(currentEdge, edgeId);
									}
									
								}
								edgeId += 1;
							}
							// std::cout << "done with multiple copies of -- " << totalProcessedEdges << " -- ";
							// printEdgeVector(currentEdge); 
							endClock = std::chrono::steady_clock::now();
							localInsertTime = std::chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
							perReportTime += localInsertTime;
						}
						else
						{
							std::cout << "Edge has less than 2 end points...\n";
						}
					}
					else
					{
						std::cout << "Edge has less than 2 end points...\n";
					}
				}
				else if(insDel == '-')
				{
					while(ss >> eEle)
					{
						tempVec.push_back(eEle);
					}
					currentEdge = tempVec;
					// if(edgeMap.find(currentEdge) != edgeMap.end())
					for(int dup = 0; dup < localAlpha; ++dup)
					{
						EdgeIdx delEdgeId = Ti[0].checkEdgeExistence(currentEdge);
						if(delEdgeId != NullEdgeIdx)
						{
							// if(totalProcessedEdges > 20700000)
							// {
							// 	std::cout << "Deleting edge " << delEdgeId << " -- "; 
							// 	printEdgeVector(currentEdge);
							// }

							edgeDeletions += 1;
							totalProcessedEdges += 1;

							beginClock = std::chrono::steady_clock::now();
							for(int copy = maxInstanceId; copy >= active; --copy)
							{
								int checkEdge = Ti[copy].checkEdgeExistenceInPendingList(currentEdge, delEdgeId);
								if(checkEdge != 0)
								{
									Ti[copy].removeEdgeFromPendingList(currentEdge, delEdgeId);
								}
								else
								{
									Ti[copy].deleteEdge(currentEdge, delEdgeId);
								}
							}

							if(active - 1 >= 1)
							{
								rho = Ti[active].getDensity();
								rhoEstActive = Ti[active].getRhoEst();

								if(rho < rhoEstActive)
								{
									active = active - 1;
									Ti[active].deleteEdge(currentEdge, delEdgeId);
								}
							}

							for(int copy = active - 1; copy >= 1; --copy)
							{
								int checkEdge = Ti[copy].checkEdgeExistenceInPendingList(currentEdge, delEdgeId);

								if(checkEdge != 0)
								{
									Ti[copy].removeEdgeFromPendingList(currentEdge, delEdgeId);
								}
								else
								{
									VertexIdx lastVertex = Ti[copy].deleteEdgeReturnLastVertex(currentEdge, delEdgeId);
									std::pair<EdgeIdx, edgeVector> eIdEdgePair = Ti[copy].getPendingEdgeForLastVertex(lastVertex);

									if(eIdEdgePair.first != NullEdgeIdx)
									{
										Ti[copy].insertEdge(eIdEdgePair.second, eIdEdgePair.first);
										Ti[copy].removeEdgeFromPendingList(eIdEdgePair.second, eIdEdgePair.first);
									}
								}

							}

							Ti[0].removeEdgeFromMap(currentEdge);

							endClock = std::chrono::steady_clock::now();
							localDeleteTime = std::chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
							perReportTime += localDeleteTime;
						}
						else
						{
							std::cout << "Edge does not exists to delete...\n";
							break;
						}
					}
					// std::cout << "done deleting multiple copies of -- " << totalProcessedEdges << " -- ";
					// printEdgeVector(currentEdge);
				}
				else if(insDel == '=')
				{
					if(Ti[active].getPendingCount() > 0)
					{
						// std::cout << "I am coming here.." << Ti[active].getPendingCount() << "\n";
						Ti[active].addListOfPendingEdges();
						// std::cout << "I am coming here.." << Ti[active].getPendingCount() << "\n";
					}

					Count yearLabel;
					ss >> yearLabel;

					std::cout << "Total Processed Edges -- " << totalProcessedEdges << "\n";
					std::cout << "Edge Additions = " << edgeAdditions << " Edge Deletions = " << edgeDeletions << "\n";

					std::vector<double> OneMinusEpsMaxInDegVector;
					std::vector<double> MaxInDegVector;
					double maxInDeg, OneMinusEpsMaxInDeg;

					MaxInDegVector.push_back(0);
					OneMinusEpsMaxInDegVector.push_back(0);

					for(int copy = 1; copy <= maxInstanceId; ++copy)
					{
						maxInDeg = (Ti[copy].getMaxLabel()*1.0)/localAlpha;
						MaxInDegVector.push_back(maxInDeg);
						OneMinusEpsMaxInDeg = Ti[copy].getDensity()/localAlpha;
						OneMinusEpsMaxInDegVector.push_back(OneMinusEpsMaxInDeg);
					}

					maxInDeg = (Ti[active].getMaxLabel()*1.0)/localAlpha;
					OneMinusEpsMaxInDeg = Ti[active].getDensity()/localAlpha;

					std::set<VertexIdx> denseSubgraph = Ti[active].querySubgraph(localAlpha);
					unsigned int denseSubgraphSize = denseSubgraph.size();
					double estimatedDensity =  Ti[active].getDensityOfInducedSubgraph(denseSubgraph)/localAlpha;

					std::cout << "Getting max partition density....\n";
					std::pair<double, unsigned int> maxPartitionedDensity = Ti[active].getMaxPartitionDensity();

					std::cout << "Max Label Val for year " << yearLabel << " : " << maxInDeg << std::endl;
					std::cout << "Max Density Val for year " << yearLabel << " : " << OneMinusEpsMaxInDeg << std::endl;
					std::cout << "Size " << yearLabel << ":" << denseSubgraphSize << "\n";
					// std::cout << "Max Indegree = " << G.getMaxIndegree() << "\n";
					std::cout << yearLabel << " Estimated Density = " << estimatedDensity << "  Size = " << denseSubgraphSize << "\n";
					std::cout << yearLabel << " Max Partition Density = " << (maxPartitionedDensity.first * 1.0)/localAlpha << "  Size = " << maxPartitionedDensity.second << "\n";
					std::cout << yearLabel << " Per Report Time = " << perReportTime << " MicroSecs \n";


					outFile << yearLabel << " " << active << " ";

					for(int copy = 1; copy <= maxInstanceId; ++copy)
					{
						outFile << copy << " " << MaxInDegVector[copy] << " " << OneMinusEpsMaxInDegVector[copy] << " ";
					}

					outFile << " " << estimatedDensity << " " << denseSubgraphSize << " " << (maxPartitionedDensity.first * 1.0)/localAlpha << " " << maxPartitionedDensity.second << " " << totalProcessedEdges << " " << edgeAdditions << " " << perReportTime <<  std::endl;

					perReportTime = 0;
				}

				// binomial sampling to scale the weight....
				// sparsify -- sampling to decide for each copy of edge to be in the graph or not...
				// maintain sparsified graph -- with edge directions....

				if(totalProcessedEdges % 10000 == 0)
				{
					std::cout << "Processed edges = " << totalProcessedEdges << "\n";
				}
			}

			// G.checkEdgeAssignment();
			// outFile.close();
		}
		else
		{
			std::cout << "First line is empty line in the file...\n";
			exit(0);
		}
		
	}
	else
	{
		std::cout << "Cannot open Graph file...\n";
		exit(0);
	}


	/*
	Graph G(v, e);

	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	G.printGraphDetails();

	// G.CountG3s();
	*/

	return 0;
}
