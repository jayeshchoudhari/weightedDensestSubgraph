#include "Graph.h"
#include <bits/stdc++.h>

using namespace std;

class Graph 
{
    private:
        Count nVertices; 								//number of vertices in the graph
        Count nEdges;     								//number of edges in this list
 		double eta = 3;
        double epsVal;

        multimap<vector<VertexIdx>, EdgeIdx> edgeMap;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        vector <vector<VertexIdx>> edgeList;

        // above could be static members shared among all the instances of the Graph...
        // each instance for a guess of density...
		
        unordered_map<EdgeIdx, VertexIdx> headOfEdgeId;

        map <VertexIdx, Count> nodeInDeg;

        vector <set<EdgeIdx>> InNbrs;							//List of InNbrs

        // ****** FOR MAINTAINING THE VISITNEXT DATA-STRUCTURE ********
        // list of neighbors for each node...
        // vector<list<VertexIdx>> listOfNeighbors;
        vector<list<EdgeIdx>> listOfNeighbors;
        // map to the elements in the neighbors list for each node...
        vector<map<EdgeIdx, list<VertexIdx>::iterator>> mapToNeighborsList;
        // vector of iterators -- one for each node... that keeps track of which element to access next in visitNext
        vector<list<EdgeIdx>::iterator> nextPositionIteratorInc;
        vector<list<EdgeIdx>::iterator> nextPositionIteratorDec;
        vector<list<EdgeIdx>::iterator> nextPositionIteratorTightInNbr;
        // *************************************************************

        // ********* TO MAINTAIN A UPDATING PRIORITY QUEUE OF OUTNEIGHBORS *************** 
        // this is to be done for each node.. and thus a vector here...
        // Each element of a vector is a map --  where the map is with a key Count, and value is a set of 
        // of edgeIds, which have the headVertices with value (indegree) as Count...
        // and nodeToOutdegMap is a reverseMap of the same...
        vector <map<Count, set<EdgeIdx>> > outdegToNodeMap;
        vector <map<EdgeIdx, Count>> nodeToOutdegMap;
        vector <unordered_map<VertexIdx, Count>> InDegreeFromNodesView;
        // ********************************************************************************

        // ********** TO MAINTAIN LABELS/NODES WITH HIGH INDEGREE *************** 
        map<Count, set<VertexIdx>> Labels;				// this is the data-structure that keeps track of top elements...
        map<VertexIdx, Count> ReverseLabels;			// this is the data-structure that keeps track of top elements...
        // **********************************************************************

    public:

		// Graph(int nv, int ne);		//gets number of vertices and edges... edges might change -- but this is just about the file... 

		Graph(Count nv, double decEta);		//gets number of vertices and edges... edges might change -- but this is just about the file... 
    	// Graph(vector <VertexIdx> &v, vector<edgeVector> &edList);		//takes in list/vector of nodes and vector of pairs which are edges (undirected)

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

		pair<VertexIdx, EdgeIdx> getTightInNbr(VertexIdx);
		VertexIdx getMinDegreeVertexInE(EdgeIdx eId);
		EdgeIdx getTightOutNbr(VertexIdx);
		EdgeIdx getMaxOutNbr(VertexIdx u);

		int updateLabels(VertexIdx u, Count changeVal);

		int insertEdge(edgeVector e, EdgeIdx eId);
		int deleteEdge(edgeVector e, EdgeIdx eId);

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

		pair<set<VertexIdx>, revItMapCountSetVertices> returnDensitySatisfiedNodes(revItMapCountSetVertices startIt, Count D);
		set<VertexIdx> getDensestSubgraph(double rgamma);
		set<VertexIdx> querySubgraph(double D_hat);
		Count getMaxIndegree();
		int checkEdgeAssignment();

		double getDensityOfInducedSubgraph(set<VertexIdx>);
		pair<double, unsigned int> getMaxPartitionDensity();

		int showPQs();
};

// Constructor to initialize the graph...
// Graph :: Graph(int nv, int ne)
Graph :: Graph(Count nv, double decEta)
{
	nVertices = nv;
	epsVal = 0.1;

	eta = decEta;

	// OutNbrs.resize(nv);

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

	// initialize adj matrix to 0...
	// each entry would be positive > 0 if an edge has some positive weight...
	/*
	for(Count i = 0; i < nv; i++)
	{
		VertexIdx v = i;
		nodeInDeg[v] = 0;
		// du[i] = 0;
		vector <Weight> eachRow;
		for(Count j = 0; j < nv; j++)
		{
			eachRow.push_back(0);
		}
		adjMatrix.push_back(eachRow);
		scaledAdjMatrix.push_back(eachRow);
		sparsedAdjMatrix.push_back(eachRow);
	}
	maxInDeg = 0;
	*/
}


int Graph :: showPQs()
{
	VertexIdx ni = 0;
	cout << "Printing outdegToNodeMap -- vector <map<Count, set<EdgeIdx>>>\n"; 
	vector <map<Count, set<EdgeIdx>>> :: iterator it1;
	for(it1 = outdegToNodeMap.begin(); it1 != outdegToNodeMap.end(); it1++)
	{
		cout << ni << "\n";
		map<Count, set<EdgeIdx>> temp = *it1;
		map<Count, set<EdgeIdx>>::iterator inIt1;
		for(inIt1 = temp.begin(); inIt1 != temp.end(); inIt1++)
		{
			cout << "---" << inIt1->first << " : "; 
			set<EdgeIdx> tset = inIt1->second;
			set<EdgeIdx>::iterator tsetIt;
			for(tsetIt = tset.begin(); tsetIt != tset.end(); tsetIt++)
			{
				cout << *tsetIt << ", ";
			}
			cout << endl;
		}
		ni++;
	}

	cout << "Printing nodeToOutdegMap -- vector <map<EdgeIdx, Count>>\n";
	ni = 0;
    vector <map<EdgeIdx, Count>>::iterator it2;
    for(it2 = nodeToOutdegMap.begin(); it2 != nodeToOutdegMap.end(); it2++)
    {
    	cout << ni << "\n";
    	map<EdgeIdx, Count> temp = *it2;
    	map<EdgeIdx, Count>::iterator inIt1;
    	for(inIt1 = temp.begin(); inIt1 != temp.end(); inIt1++)
    	{
    		cout << " --- " << inIt1->first << " " << inIt1->second << "\n"; 
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

Count Graph :: getNumVertices()
{
	return nVertices;
}

/*
Count Graph :: getNumEdges()
{
	return nEdges;
}
*/

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

	// cout << nodeInDeg[w] << " headNode indegree during addition\n";
	addDirectedEdgeToInOutNbrs(eId, w);

	// cout << InNbrs[w].size() << " headNode innbrs size during addition\n";

	// headOfEdgeId[eId] = w;

	// check if this results into making some neighboring edge of w tight...
	pair<VertexIdx, EdgeIdx> minNbrNodeEdgePair = getTightInNbr(w);
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

	VertexIdx headNode = headOfEdgeId[eId];

	cout << nodeInDeg[headNode] << " headNode indegree before deletion\n";
	cout << InNbrs[headNode].size() << " headNode innbrs size before deletion\n";
	
	// remove e from the InNbrs of headNode...
	removeDirectedEdgeFromInOutNbrs(eId, headNode);
	cout << nodeInDeg[headNode] << "after deletion headNode indegree before deletion -- shouldnt change yet..\n";

	// check if u belongs to in-neighbors of v
	
	/*
	if(InNbrs[v][u] != 0)
	{
		eTupleUnWeighted eDelete(u,v);
		removeDirectedEdge(eDelete);
		w = v;
	}
	else
	{
		eTupleUnWeighted eDelete(v,u);
		removeDirectedEdge(eDelete);	
		w = u;
	}
	*/
	lastEId = eId;
	w = headNode;
	ePrime = getTightOutNbr(w);
	
	map<EdgeIdx, int> flippedEdges;

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

	return 0;
}


int Graph :: updateIncPointer(VertexIdx u)
{
	// shift the itertor to the next position...
	list<VertexIdx>:: iterator updateItInc;

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
	list<VertexIdx>:: iterator updateItDec;

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
	list<VertexIdx>:: iterator updateItNext;

	updateItNext = nextPositionIteratorTightInNbr[u];
	updateItNext++;
	if(updateItNext == listOfNeighbors[u].end())
	{
		updateItNext = listOfNeighbors[u].begin();
	}
	nextPositionIteratorTightInNbr[u] = updateItNext;

	return 0;
}




// int Graph :: addToInNbrs(VertexIdx u, VertexIdx v)
int Graph :: addEdgeToInNbrsForVisitNext(VertexIdx headNode, EdgeIdx eId)
{
	// This is to add v to the in-neighbors of u...
	// add to the list of in-neighbors and update the existence and 
	// address of this new neigbor in the map...
	listOfNeighbors[headNode].push_back(eId);				// adding v to the in-neighbors of u...
	list<EdgeIdx>::iterator esAddress = listOfNeighbors[headNode].end();
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


// int Graph :: addDirectedEdge(eTupleUnWeighted e)
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


// int Graph :: flipDirectedEdge(eTupleUnWeighted e)
int Graph :: flipDirectedEdge(EdgeIdx eId, VertexIdx oldHeadNode, VertexIdx newHeadNode)
{
	removeDirectedEdgeFromInOutNbrs(eId, oldHeadNode);

	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	// eTupleUnWeighted flippedEdge (v,u);
	addDirectedEdgeToInOutNbrs(eId, newHeadNode);
	
	return 0;
}

int Graph :: updateNextNeighbors(VertexIdx headNode, Count newDuVal, int incOrDec)
{
	list<EdgeIdx> :: iterator updateItInc;
	unordered_map<EdgeIdx, int> touchedNeighbors;

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


pair<VertexIdx, EdgeIdx> Graph :: getTightInNbr(VertexIdx u)
{
	// VertexIdx neighborToReturn = NullVertexIdx;
	EdgeIdx neighborEdgeIdToReturn = NullEdgeIdx;
	VertexIdx neighborNodeIdToReturn = NullVertexIdx;

	list<EdgeIdx> :: iterator updateItInc;

	Count start = 0;
	Count maxNumNeighborsToCheck = (Count)(4 * nodeInDeg[u]) / eta;
	
	unordered_map<EdgeIdx, int> touchedNeighbors;
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

	return make_pair(neighborNodeIdToReturn, neighborEdgeIdToReturn);
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
	// get the degToNode map of u
	if(outdegToNodeMap[u].size() >= 1)
	{
		map<Count, set<EdgeIdx>>::reverse_iterator rit = outdegToNodeMap[u].rbegin();
		Count maxVal = rit->first;
		set<EdgeIdx> maxValSet = rit->second;
		set<EdgeIdx>::iterator it = maxValSet.begin();
		EdgeIdx maxEle = *it;

		return maxEle;
	}

	return NullVertexIdx;
}


EdgeIdx Graph :: checkEdgeExistence(edgeVector e)
{
	multimap<vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
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
	edgeMap.insert(pair<vector<VertexIdx>, EdgeIdx>(e, eId));
	return 0;
}

EdgeIdx Graph :: getEdgeId(edgeVector e)
{
	multimap<vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
	EdgeIdx eId = edgeMapIt->second;
	return eId;
}

int Graph :: removeEdgeFromMap(edgeVector e)
{
	// edgeMap.erase(e);
	multimap<vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt = edgeMap.find(e);
	edgeMap.erase(edgeMapIt);
	return 0;
}


double Graph :: getDensity()
{
	map<Count, set<VertexIdx>>::reverse_iterator labelsIt = Labels.rbegin();
	double currMaxDensityValue = (labelsIt->first)*(1-epsVal);
	return currMaxDensityValue;
}

Count Graph :: getLabel(VertexIdx u)
{
	return ReverseLabels[u];
}

Count Graph :: getMaxLabel()
{
	map<Count, set<VertexIdx>> :: reverse_iterator rit = Labels.rbegin();
	Count maxVal = rit->first;
	/*
	set<VertexIdx> maxValSet = rit->second;
	set<VertexIdx>::iterator it = maxValSet.begin();
	VertexIdx maxEle = *it;
	*/
	return maxVal;
}

pair<set<VertexIdx>, revItMapCountSetVertices> Graph :: returnDensitySatisfiedNodes(revItMapCountSetVertices startIt, Count D)
{
	set<VertexIdx> B;
    map<Count, set<VertexIdx>>::reverse_iterator preservedRit;
    map<Count, set<VertexIdx>>::reverse_iterator rit;

	for(rit = startIt; rit != Labels.rend(); ++rit)
	{
		preservedRit = rit;
		Count densityVal = rit->first;
		if(densityVal >= D)
		{
			set<VertexIdx> dvertices = rit->second;
			set<VertexIdx> :: iterator dvIt;
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

	pair rP = make_pair(B, preservedRit);
    return rP;
}


set<VertexIdx> Graph :: getDensestSubgraph(double rgamma)
{
	Count decrementVal = 1;

	double D = getMaxLabel();
	set<VertexIdx> A, B;
	unsigned int ASize, BSize;
	double sizeRatio;

	revItMapCountSetVertices rit;

	rit = Labels.rbegin();
	pair<set<VertexIdx>, revItMapCountSetVertices> AElementsPair = returnDensitySatisfiedNodes(rit, D);
	A = AElementsPair.first;

	if(AElementsPair.second != Labels.rend())
	{
		rit = AElementsPair.second;
	}
	else
	{
		cout << "Returning just the elements in A -- which is a set of all elements...\n";
		return A;
	}

	// D = D-eta;
	D = D-decrementVal;
	B = A;

	pair<set<VertexIdx>, revItMapCountSetVertices> newElementsToBPair = returnDensitySatisfiedNodes(rit, D);
	set<VertexIdx> newElementsToB = newElementsToBPair.first;

	rit = newElementsToBPair.second;

	set<VertexIdx> :: iterator bIt;	

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
			cout << "Cannot accumulate B anymore...\n";
			return B;
		}
	}

	return B;
}


set<VertexIdx> Graph :: querySubgraph(double D_hat)
{
	double rgamma = sqrt(2 * eta * log(nVertices) / D_hat);
	cout << "rGamma = " << rgamma << endl;
	return getDensestSubgraph(rgamma);
}

Count Graph :: getMaxIndegree()
{
	map<VertexIdx, Count> :: iterator it;
	Count maxVal = 0;
	VertexIdx maxNode;
	for(it = nodeInDeg.begin(); it != nodeInDeg.end(); ++it)
	{
		if(it->second > maxVal)
		{
			maxVal = it->second;
			maxNode = it->first;
		}
	}
	cout << "maxNode -- " << maxNode << "----"; 
	return maxVal;
}



int Graph :: checkEdgeAssignment()
{

	cout << "Checking Edge Assignment --- " << endl;
	Count i = 0;
	multimap<vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt;

	for(edgeMapIt = edgeMap.begin(); edgeMapIt != edgeMap.end(); ++edgeMapIt)
	{
		i++;
		vector<VertexIdx> evector = edgeMapIt->first;
		EdgeIdx eId = edgeMapIt->second;

		VertexIdx headNode = headOfEdgeId[eId];
		Count headInDeg = nodeInDeg[headNode];

		VertexIdx minDegNode = getMinDegreeVertexInE(eId);
		Count minInDeg = nodeInDeg[minDegNode];

		Count diffVal = headInDeg - minInDeg;

		if(diffVal >= 4)
		{
			cout << eId << " " << diffVal << "\n"; 
		} 		
	}

	cout << "Checked for " << i << " edges --- " << endl;

	return 0;
}

double Graph :: getDensityOfInducedSubgraph(set<VertexIdx> denseSubgraphNodes)
{

	cout << "Checking density of said nodes... -- ";
	Count i = 0;
	multimap<vector<VertexIdx>, EdgeIdx>::iterator edgeMapIt;

	vector<edgeVector> selectedEdges;

	Count edgeCount = 0;

	vector<VertexIdx> evector;

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

	cout << i << "---" << (edgeCount*1.0)/denseSubgraphNodes.size() << " " << denseSubgraphNodes.size() << endl;

	return estDensity;
}


pair<double, unsigned int> Graph :: getMaxPartitionDensity()
{
	double maxDensity = 0;
	unsigned int maxSubgraphSize = 0;

	set<VertexIdx> denseSubgraphNodes;

	revItMapCountSetVertices rit;

	for(rit = Labels.rbegin(); rit != Labels.rend(); ++rit)
	{
		set<VertexIdx> newElementsToSet = rit->second;

		set<VertexIdx> :: iterator bIt;
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

	return make_pair(maxDensity, maxSubgraphSize);
}


// utilities...
int sampleFromBinomial(int wt, double p)
{
	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// default_random_engine generator (seed);
  	random_device rd;
    mt19937 gen(rd());
    // perform 4 trials, each succeeds 1 in 2 times
    // std::binomial_distribution<> d(4, 0.5);

  	std::binomial_distribution<int> distribution (wt, p);
  	return distribution(gen);
}


int printEdgeVector(edgeVector e)
{
	for(unsigned int i = 0; i < e.size(); i++)
	{
		cout << e[i] << " ";
	}
	cout << "\n";

	return 0;
}

int main(int argc, char** argv)
{
	// std::vector <VertexIdx> v = {0,1,2,3};
	// std::vector<pair<VertexIdx, VertexIdx>> e = {{1, 3}, {3, 0}, {0, 2}, {2, 1}, {2, 3}};

	// each file starts with n, m, and max_weight i.e. two numbers on each line
	// after that there are m lines each with the edge information
	// each line is a 4 tuple -- src dest weight 1/0
	// last element of tuple: 1 - insert, 0 - delete 
	// if last element of tuple is 0, check if it was inserted before -- checkEdgeExistence()
	// after inserting each edge the graph is updated...

	//read each line from file...
	// string graphFileName = "dblp.theory.hg.dyn.txt";
	// string graphFileName = "dblp.theory.hg.dyn.10.txt";
	// string graphFileName = "dblp.theory.hg.dyn.15.txt";
	// string graphFileName = "dblp.theory.hypergraph.txt";
	// string graphFileName = "sampleGraph-1.txt";
	// string graphFileName = "sample-hypergraph-1.txt";

	string graphFileName = argv[1];
	string outFileName = argv[2];


	ifstream graphFile;
	ofstream outFile;

	string line;
	stringstream ss;
	
	Count n, m, maxKEdge, totalProcessedEdges = 0; 
	VertexIdx e_src, e_dest, eEle;
	Weight e_weight, Wmax;
	char insDel;

	double scalingProbParam_c, sparsifyProbParam_c, scalingEps, sparsifyEps;
	double scalingEpsSquared = scalingEps * scalingEps;
	double sparsifyEpsSquared = sparsifyEps * sparsifyEps;

	double scalingProb, minWeightedDensity; 

	// cout << "Did I come here...\n";

	graphFile.open(graphFileName, ifstream::in);
	if(graphFile.is_open())
	{
		// if (getline(graphFile, line) != -1)
		if (getline(graphFile, line))
		{
			outFile.open(outFileName, ofstream::out);

			outFile << "YearLabel" << " " << "MaxInDeg" << " " << "(1-eps)MaxInDeg" << " " << "EstimatedDensity" << " " << "DenseSubgraphSize" << " " << "MaxPartitionDensity" << " " << "MaxPartitionSize" << endl;

			ss.clear();
			ss.str("");
			ss << line;
			// ss >> n >> m >> Wmax;
			ss >> n >> maxKEdge;

			// for now, let us say that only the exact number of nodes is only known... 
			// initialize the datastructure required for each node...
			// d(u) 
			// InNbrs -- List of in-nbrs for each node
			// OutNbrs -- Max priority queue for each node u... indexed by d(v)

			double etaVal = stod(argv[3]);

			Graph G(n, etaVal); 	// This initializes du, InNbrs, OutNbrs(list of priority queues...)
			cout << "Done with initialization...\n";
		

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
							edgeAdditions += 1;
							totalProcessedEdges += 1;
							G.setEdgeId(currentEdge, edgeId);
							G.addEdgeToEdgeList(currentEdge);
							
							cout << "Add the edge -- " << edgeId << " -- "; 
							printEdgeVector(currentEdge);
							G.insertEdge(currentEdge, edgeId);

							edgeId += 1;
						}
						else
						{
							cout << "Edge has less than 2 end points...\n";
						}
					}
					else
					{
						cout << "Edge has less than 2 end points...\n";
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
					EdgeIdx delEdgeId = G.checkEdgeExistence(currentEdge);
					if(delEdgeId != NullEdgeIdx)
					{
						edgeDeletions += 1;
						totalProcessedEdges += 1;
						cout << "Deleting edge " << delEdgeId << " -- "; 
						printEdgeVector(currentEdge);
						G.deleteEdge(currentEdge, delEdgeId);
						G.removeEdgeFromMap(currentEdge);
					}
					else
					{
						cout << "Edge does not exists to delete...\n";
					}
				}
				else if(insDel == '=')
				{
					Count yearLabel;
					ss >> yearLabel;

					cout << "Total Processed Edges -- " << totalProcessedEdges << "\n";
					cout << "Edge Additions = " << edgeAdditions << " Edge Deletions = " << edgeDeletions << "\n";

					double maxInDeg = G.getMaxLabel();
					double OneMinusEpsMaxInDeg = G.getDensity();
					set<VertexIdx> denseSubgraph = G.querySubgraph(1);
					unsigned int denseSubgraphSize = denseSubgraph.size();
					double estimatedDensity =  G.getDensityOfInducedSubgraph(denseSubgraph);

					cout << "Getting max partition density....\n";
					pair<double, unsigned int> maxPartitionedDensity = G.getMaxPartitionDensity();

					cout << "Max Label Val for year " << yearLabel << " : " << maxInDeg << endl;
					cout << "Max Density Val for year " << yearLabel << " : " << OneMinusEpsMaxInDeg << endl;
					cout << "Size " << yearLabel << ":" << denseSubgraphSize << "\n";
					// cout << "Max Indegree = " << G.getMaxIndegree() << "\n";
					cout << yearLabel << " Estimated Density = " << estimatedDensity << "  Size = " << denseSubgraphSize << "\n";
					cout << yearLabel << " Max Partition Density = " << maxPartitionedDensity.first << "  Size = " << maxPartitionedDensity.second << "\n";

					outFile << yearLabel << " " << maxInDeg << " " << OneMinusEpsMaxInDeg << " " << estimatedDensity << " " << denseSubgraphSize << " " << maxPartitionedDensity.first << " " << maxPartitionedDensity.second << endl;
				}

				// binomial sampling to scale the weight....
				// sparsify -- sampling to decide for each copy of edge to be in the graph or not...
				// maintain sparsified graph -- with edge directions....

				/*
				Weight scaledWeight = sampleFromBinomial(e_weight, scalingProb);

				for(Weight w = 0; w < scaledWeight; w++)
				{

				}
				*/
			}

			G.checkEdgeAssignment();
			outFile.close();

		}
		else
		{
			cout << "First line is empty line in the file...\n";
			exit(0);
		}
		
	}
	else
	{
		cout << "Cannot open Graph file...\n";
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