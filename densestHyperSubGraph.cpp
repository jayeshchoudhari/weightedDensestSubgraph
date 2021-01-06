#include "Graph.h"
#include <bits/stdc++.h>

using namespace std;
// using namespace rwNamespace;

class Graph 
{
    private:
        Count nVertices; 								//number of vertices in the graph
        Count nEdges;     								//number of edges in this list
        VertexIdx NullVertexIdx = -1;
        EdgeIdx NullEdgeIdx = -1;

        vector <VertexIdx> srcs;      					//array of source vertices
        vector <VertexIdx> dsts;      					//array of destination vertices
        
        vector <VertexIdx> nodeList;
        vector <vector<VertexIdx>> adjList;				//adj List
        vector <vector<Weight>> adjMatrix;				//adj Matrix
        vector <vector<Weight>> scaledAdjMatrix;		//scaled adj Matrix
        vector <vector<Weight>> sparsedAdjMatrix;		//sparsed adj Matrix

        unordered_map<vector<VertexIdx>, Count> edgeMap;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        unordered_map<Count, VertexIdx> headOfEdgeId;

        vector <vector<VertexIdx>> edgeList;
        unordered_map <VertexIdx, Count> nodeInDeg;

        vector <Count> du;
        // vector <vector<Count>> InNbrs;					//List of InNbrs
        vector <set<EdgeIdx>> InNbrs;					//List of InNbrs
        vector <priority_queue<Count> > OutNbrs;
        vector <Count> nextNeighbor;

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
        // *************************************************************


        // ********** TO MAINTAIN LABELS/NODES WITH HIGH INDEGREE *************** 
        map<Count, set<VertexIdx>> Labels;				// this is the data-structure that keeps track of top elements...
        map<VertexIdx, Count> ReverseLabels;			// this is the data-structure that keeps track of top elements...
        // *************************************************************

        Count maxLabel;					//Node with the max Indegree


        float eta;

    public:

		// Graph(int nv, int ne);		//gets number of vertices and edges... edges might change -- but this is just about the file... 

		Graph(int nv);		//gets number of vertices and edges... edges might change -- but this is just about the file... 

    	Graph(vector <VertexIdx> &v, vector<eTupleUnWeighted> &edList);		//takes in list/vector of nodes and vector of pairs which are edges (undirected)

    	int printGraphDetails();

  		// int CountG3s();
		// eTupleUnWeighted lStepRandomWalk(Count, VertexIdx);

		Count getNumVertices();
		Count getNumEdges();
        // //Make a copy
        // Graph copy() const;

        // //For debugging.  Not for serialization.
        // void print(FILE *f = stdout) const;

		// int initializeDu(Count n);

		// int addDirectedEdge(eTupleUnWeighted);
		int addDirectedEdgeToInOutNbrs(EdgeIdx eId, VertexIdx v);
		int removeDirectedEdgeFromInOutNbrs(EdgeIdx eId);
		int flipDirectedEdge(eTupleUnWeighted);
		// int addToPriorityQueue(VertexIdx, VertexIdx, Count);
		int addToPriorityQueue(VertexIdx, EdgeIdx, Count);
		int removeFromPriorityQueue(VertexIdx, EdgeIdx);
		int updateNextNeighbors(VertexIdx u, Count newDuVal, int incOrDec);

		int incrementDu(VertexIdx);
		int decrementDu(VertexIdx);

		EdgeIdx getTightInNbr(VertexIdx);
		VertexIdx getTightOutNbr(VertexIdx);

		Count getLabel(VertexIdx);
		Count getMaxLabel();
		int updateLabels(VertexIdx u, Count changeVal);
		VertexIdx getMaxOutNbr(VertexIdx u);

		int insertEdge(edgeVector, int);
		int deleteEdge(edgeVector, int);

		int updateIncPointer(VertexIdx u);
		int updateDecPointer(VertexIdx u);
		int updateTightInNbrIterator(VertexIdx u);
		int addToInNbrs(VertexIdx u, VertexIdx v);
		int removeFromInNbrs(VertexIdx u, VertexIdx v);
};

// Constructor to initialize the graph...
// Graph :: Graph(int nv, int ne)
Graph :: Graph(int nv)
{
	nVertices = nv;
	OutNbrs.resize(nv);
	du.resize(nv);
	outdegToNodeMap.resize(nv);
	nodeToOutdegMap.resize(nv);
	// nextNeighbor.resize(nv);

	// to maintain visitNext structure... 
	listOfNeighbors.resize(nv);
	mapToNeighborsList.resize(nv);
	nextPositionIteratorInc.resize(nv);
	nextPositionIteratorDec.resize(nv);
	nextPositionIteratorTightInNbr.resize(nv);

	// initialize adj matrix to 0...
	// each entry would be positive > 0 if an edge has some positive weight...
	for(Count i = 0; i < nv; i++)
	{
		du[i] = 0;
		vector <Weight> eachRow;
		for(Count j = 0; j < nv; j++)
		{
			eachRow.push_back(0);
		}
		adjMatrix.push_back(eachRow);
		scaledAdjMatrix.push_back(eachRow);
		sparsedAdjMatrix.push_back(eachRow);
		// InNbrs.push_back(eachRow);
		// nextNeighbor[i] = 0;
	}
}


Graph :: Graph(vector <VertexIdx> &v, vector<eTupleUnWeighted> &e)
// Graph::Graph(vector <VertexIdx> &v, vector<eTupleUnWeighted> &edList)
{
	nVertices = v.size();
	nEdges = e.size();

	// std::cout << "num of nodes = " << nVertices << " num of edges = " << nEdges << "\n";

	vector<VertexIdx>::iterator vIt;

	for(vIt = v.begin(); vIt != v.end(); vIt++)
	{
		nodeList.push_back(*vIt);
	}
	
	// std::cout << "got vertices \n"; 

	adjList.resize(nVertices);

	vector<eTupleUnWeighted>::iterator edgeIt;

	// for(const eTupleUnWeighted &edgeIt : e)
	for(edgeIt = e.begin(); edgeIt != e.end(); edgeIt++)
	{
		VertexIdx src = get<0>(*edgeIt);
		VertexIdx dest = get<1>(*edgeIt);

		// populate adjacency list
		adjList[src].push_back(dest);
		adjList[dest].push_back(src);

		// populate list of edges
		vector <VertexIdx> tempEdge;
		tempEdge.push_back(src);
		tempEdge.push_back(dest);
		edgeList.push_back(tempEdge);

		// populate node degress
		nodeInDeg[src] += 1;
		nodeInDeg[dest] += 1;
	}
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

Count Graph :: getNumEdges()
{
	return nEdges;
}

int Graph :: printGraphDetails()
{
	std::cout << "num of nodes = " << nVertices << " num of edges = " << nEdges << "\n";

	unordered_map<VertexIdx, Count>::iterator nIt;
	for(nIt = nodeInDeg.begin(); nIt != nodeInDeg.end(); nIt++)
	{
		std::cout << nIt->first << " : " << nIt->second << "\n";
	}

	
	// Print graph  //
	for (int i = 0; i < nVertices; i++)
	{
		// print current vertex number
		cout << i << " --> ";

		// print all neighboring vertices of vertex i
		for (int v : adjList[i])
			cout << v << " ";
		cout << endl;
	}
	// 
	return 0;
}

/*
int Graph :: checkEdgeExistence(e)
{
	// check in the adjacency list... if the pair exists....
	return flag;
}
*/

/*
int Graph :: initializeDu(Count numVertices)
{
	for(VertexIdx i = 0; i < numVertices; i++)
	{
		du[i] = 0;
	}
	return 0;
}
*/

int Graph :: addToPriorityQueue(VertexIdx u, EdgeIdx eId, Count vVal)
{
	// we need to add/update v in the priority queue of u;
	// Count oldVal = nodeToOutdegMap[u][eId];
	nodeToOutdegMap[u][eId] = vVal;

	// outdegToNodeMap[u][oldVal].erase(eId);
	outdegToNodeMap[u][vVal].insert(eId);

	return 0;
}


int Graph :: removeFromPriorityQueue(VertexIdx u, EdgeIdx eId)
{
	// remove/decrement v in the priority queue of u;
	Count oldVal = nodeToOutdegMap[u][eId];			// gets the value of the degree of head node of eId...
	// Count newVal = oldVal - 1;
	nodeToOutdegMap[u].erase(eId);
	outdegToNodeMap[u][oldVal].erase(eId);
	// outdegToNodeMap[u][newVal].insert(eId);

	return 0;
}


int Graph :: insertEdge(edgeVector e, int eId)
{
	// e = u,v -- at this moment an edge isnt directed...
	// hereafter it will be... from here onwards we would orient the edges...
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	VertexIdx w, wPrime;
	EdgeIdx ePrime;

	edgeVector::iterator minIt = min_element(e.begin(), e.end());
	VertexIdx minDegVertex = distance(e.begin(), minIt);

	w = minDegVertex;

	addDirectedEdgeToInOutNbrs(eId, w);

	// headOfEdgeId[eId] = w;

	// check if this results into making some neighboring edge of w tight... 
	ePrime = getTightInNbr(w);		// tight edge... -- head of this edge give tight neighbor node...
	// while(wPrime != -1)
	while(ePrime != -1)
	{
		wPrime = headOfEdgeId[ePrime];
		// eTupleUnWeighted eFlip(wPrime, w);
		flipDirectedEdge(ePrime, w);
		w = wPrime;
		wPrime = getTightInNbr(w);
	}

	incrementDu(w);

	incrementDu(w);
	return 0;
}


int Graph :: deleteEdge(eTupleUnWeighted e)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	VertexIdx w, wPrime;
	edgeVector::iterator minIt = min_element(e.begin(), e.end());
	VertexIdx minDegVertex = distance(e.begin(), minIt);

	// check if u belongs to in-neighbors of v
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

	
	while(getTightOutNbr(w) != NullVertexIdx)
	{
		wPrime = getTightOutNbr(w);
		eTupleUnWeighted eFlip(w, wPrime);
		flipDirectedEdge(eFlip);
		wPrime = w;
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
int Graph :: addToInNbrs(VertexIdx headNode, EdgeIdx eId)
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
	}

	return 0;
}


// int Graph :: removeFromInNbrs(VertexIdx u, VertexIdx v)
int Graph :: removeFromInNbrs(VertexIdx headNode, EdgeIdx eId)
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
			updateIncPointer(headNode, eId);
		}

		if(nextPositionIteratorDec[headNode] == mapToNeighborsList[headNode][eId])
		{
			updateDecPointer(headNode, eId);
		}

		if(nextPositionIteratorTightInNbr[headNode] == mapToNeighborsList[headNode][eId])
		{
			updateTightInNbrIterator(headNode, eId);
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
int Graph :: addDirectedEdgeToInOutNbrs(EdgeIdx eId, VertexIdx v);
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	VertexIdx headNode = v;
	headOfEdgeId[eId] = headNode;

	// add u to in-neighbors of v
	// InNbrs[v][u] += 1;
	InNbrs[headNode].insert(eId);
	addToInNbrs(headNode, eId);

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
	Count dvVal = nodeInDeg[headNode];

	// update e in the priority queue of all u's except v with this v's value...
	edgeVector e = edgeList[eId];
	for(unsigned int i = 0; i < e.size(); i++)
	{
		VertexIdx u = e[i];
		if(u != headNode)
		{
			addToPriorityQueue(u, eId, dvVal);
		}
	}

	return 0;
}


int Graph :: removeDirectedEdgeFromInOutNbrs(EdgeIdx eId)
{
	// e = u,v directed
	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);

	VertexIdx headNode = headOfEdgeId[eId];

	// remove/decrement u from in-neighbors of v
	// InNbrs[v][u] -= 1;

	InNbrs[headNode].erase(eId);
	removeFromInNbrs(eId, headNode);

	// Note that this is to be done only if the u no more remains an in-neighbor of v..
	// remember this is a multigraph.... 
	// if(InNbrs[v][u] == 0)
	// {
	// 	removeFromInNbrs(v, u);
	// }

	// remove/decrement v from the priority queue out-neighbors of u
	edgeVector e = edgeList[eId];
	for(unsigned int i = 0; i < e.size(); i++)
	{
		VertexIdx u = e[i];
		if(u != headNode)
		{
			removeFromPriorityQueue(u, eId);
		}
	}
	return 0;
}


// int Graph :: flipDirectedEdge(eTupleUnWeighted e)
int Graph :: flipDirectedEdge(edgeVector e, VertexIdx changeHeadToNode)
{
	removeDirectedEdgeFromInOutNbrs(e);

	// VertexIdx u = get<0>(e);
	// VertexIdx v = get<1>(e);
	// eTupleUnWeighted flippedEdge (v,u);
	addDirectedEdgeToInOutNbrs(e, changeHeadToNode);
	
	return 0;
}

int Graph :: updateNextNeighbors(VertexIdx u, Count newDuVal, int incOrDec)
{
	list<VertexIdx> :: iterator updateItInc;
	unordered_map<VertexIdx, int> touchedNeighbors;

	Count start = 0;
	Count maxNumNeighborsToUpdate = (Count) (4 * nodeInDeg[u]) / eta;
	// note that you shouldnt be updating an element multiple times within a same 
	// update... so keep track of the updated neighbors... 
	// every time start with a new/empty list of neighbors that would be updated...
	if(incOrDec == 1)
	{
		updateItInc = nextPositionIteratorInc[u];
	}
	else if (incOrDec == -1)
	{
		updateItInc = nextPositionIteratorDec[u];
	}

	while(start < maxNumNeighborsToUpdate)
	{
		// access the next neighbor 
		VertexIdx usNextNeighbor = *updateItInc;
		
		if(touchedNeighbors.find(usNextNeighbor) == touchedNeighbors.end())
		{		
			// update new in-degree value of u to the neighbor
			addToPriorityQueue(usNextNeighbor, u, newDuVal);
	
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
		nextPositionIteratorInc[u] = updateItInc;
	}
	else
	{
		nextPositionIteratorDec[u] = updateItInc;
	}

	return 0;
}

int Graph :: incrementDu(VertexIdx u)
{
	updateLabels(u, 1);
	nodeInDeg[u] += 1;
	Count newDuVal = nodeInDeg[u];

	// update 4 din(u)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(u, newDuVal, 1);

	return 0;
}

int Graph :: decrementDu(VertexIdx u)
{
	updateLabels(u, -1);
	nodeInDeg[u] -= 1;
	Count newDuVal = nodeInDeg[u];

	// update 4 din(u)/eta next in-neigbors of u about the change in the in-degree of u 
	updateNextNeighbors(u, newDuVal, -1);
	
	return 0;
}

int Graph :: updateLabels(VertexIdx u, Count changeVal)
{
	// updating the Labels data-structure...
	Labels[nodeInDeg[u]].erase(u);
	nodeInDeg[u] +=  changeVal;
	Labels[nodeInDeg[u]].insert(u);
	ReverseLabels[u] = nodeInDeg[u];

	return 0;
}

EdgeIdx Graph :: getTightInNbr(VertexIdx u)
{
	// VertexIdx neighborToReturn = NullVertexIdx;
	EdgeIdx neighborEdgeIdToReturn = NullEdgeIdx;

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
			VertexIdx edgeHeadNode = headOfEdgeId[usNextNeighborEdgeId]; 
			if(nodeInDeg[edgeHeadNode] <= nodeInDeg[u] - eta/2)
			{
				neighborEdgeIdToReturn = usNextNeighborEdgeId;
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

	return neighborEdgeIdToReturn;
}


VertexIdx Graph :: getTightOutNbr(VertexIdx u)
{
	VertexIdx t = getMaxOutNbr(u);
	if(nodeToOutdegMap[u][t] >= nodeInDeg[u] + eta/2)
	{
		return t;
	}
	return NullVertexIdx;
}


Count Graph :: getLabel(VertexIdx u)
{
	return ReverseLabels[u];
}

Count Graph :: getMaxLabel()
{
	map<Count, set<VertexIdx>> :: reverse_iterator rit = Labels.rbegin();
	Count maxVal = rit->first;

	set<VertexIdx> maxValSet = rit->second;
	set<VertexIdx>::iterator it = maxValSet.begin();
	VertexIdx maxEle = *it;

	return maxVal;
}


VertexIdx Graph :: getMaxOutNbr(VertexIdx u)
{
	// get the degToNode map of u
	map<Count, set<VertexIdx>>::reverse_iterator rit = outdegToNodeMap[u].rbegin();
	Count maxVal = rit->first;
	set<VertexIdx> maxValSet = rit->second;
	set<VertexIdx>::iterator it = maxValSet.begin();
	VertexIdx maxEle = *it;

	return maxEle;
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



int main()
{
	std::vector <VertexIdx> v = {0,1,2,3};
	std::vector<pair<VertexIdx, VertexIdx>> e = {{1, 3}, {3, 0}, {0, 2}, {2, 1}, {2, 3}};

	// each file starts with n, m, and max_weight i.e. two numbers on each line
	// after that there are m lines each with the edge information
	// each line is a 4 tuple -- src dest weight 1/0
	// last element of tuple: 1 - insert, 0 - delete 
	// if last element of tuple is 0, check if it was inserted before -- checkEdgeExistence()
	// after inserting each edge the graph is updated...

	//read each line from file...
	string graphFileName = "dblp.theory.hypergraph.txt";
	ifstream graphFile;

	string line;
	stringstream ss;
	
	Count n, m, maxKEdge;
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

			Graph G(n); 	// This initializes du, InNbrs, OutNbrs(list of priority queues...)
			cout << "Done with initialization...\n";
			

			minWeightedDensity = (Wmax*1.0)/2;

			// scalingProb = (c log n) / (eps^2 * rho_min)  
			scalingProb = (scalingProbParam_c * log2(n)) / (scalingEpsSquared * minWeightedDensity);

			Count edgeId = 0;

			while(getline(graphFile, line))
			{
				ss.clear();
				ss.str("");
				ss << line;

				ss >> insDel;

				hyperEdge tempVec;

				while(ss >> eEle)
				{
					tempVec.push_back(eEle);
				}

				//pop_back
				Count yearVal = tempVec[tempVec.size()-1]; 	// last element is the year/timestamp value...
				tempVec.pop_back();

				edgeVector currentEdge = tempVec;

				if(insDel == "+")
				{
					if(edgeMap.find(currentEdge) == edgeMap.end())
					{
						edgeMap[currentEdge] = edgeId;
						G.addEdgeToEdgeList(currentEdge);
						cout << "Add the edge -- " << edgeId << endl;

						insertEdge(currentEdge, edgeId);

						edgeId += 1;
					}
					else
					{
						cout << "Edge already exists in the graph...\n";
					}
				}
				else if(insDel == "-")
				{
					if(edgeMap.find(currentEdge) == edgeMap.end())
					{
						cout << "Deleting edge " << edgeMap[currentEdge] << endl;
					}
					else
					{
						cout << "Edge does not exists to delete...\n";
					}
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