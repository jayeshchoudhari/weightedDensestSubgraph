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

        vector <VertexIdx> srcs;      					//array of source vertices
        vector <VertexIdx> dsts;      					//array of destination vertices
        
        vector <VertexIdx> nodeList;
        vector <vector<VertexIdx>> adjList;				//adj List
        vector <vector<Weight>> adjMatrix;				//adj Matrix
        vector <vector<Weight>> scaledAdjMatrix;		//scaled adj Matrix
        vector <vector<Weight>> sparsedAdjMatrix;		//sparsed adj Matrix

        vector <vector<VertexIdx>> edgeList;
        unordered_map <VertexIdx, Count> nodeInDeg;

        vector <Count> du;
        vector <vector<Count>> InNbrs;					//List of InNbrs
        vector <priority_queue<Count> > OutNbrs;
        vector <Count> nextNeighbor;

        // ****** FOR MAINTAINING THE VISITNEXT DATA-STRUCTURE ********
        // list of neighbors for each node...
        vector<list<VertexIdx>> listOfNeighbors;
        // map to the elements in the neighbors list for each node...
        vector<map<VertexIdx, list<VertexIdx>::iterator>> mapToNeighborsList;
        // vector of iterators -- one for each node... that keeps track of which element to access next in visitNext
        vector<list<VertexIdx>::iterator> nextPositionIteratorInc;
        vector<list<VertexIdx>::iterator> nextPositionIteratorDec;
        vector<list<VertexIdx>::iterator> nextPositionIteratorTightInNbr;
        // *************************************************************

        // ********* TO MAINTAIN A UPDATING PRIORITY QUEUE OF OUTNEIGHBORS *************** 
        // this is to be done for each node.. and thus a vector here...
        vector <map<Count, set<VertexIdx>> > outdegToNodeMap;
        vector <map<VertexIdx, Count>> nodeToOutdegMap;
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

		int addDirectedEdge(eTupleUnWeighted);
		int removeDirectedEdge(eTupleUnWeighted);
		int flipDirectedEdge(eTupleUnWeighted);
		int addToPriorityQueue(VertexIdx, VertexIdx, Count);
		int removeFromPriorityQueue(VertexIdx, VertexIdx);
		int updateNextNeighbors(VertexIdx u, Count newDuVal, int incOrDec);

		int incrementDu(VertexIdx);
		int decrementDu(VertexIdx);

		VertexIdx getTightInNbr(VertexIdx);
		VertexIdx getTightOutNbr(VertexIdx);

		Count getLabel(VertexIdx);
		Count getMaxLabel();
		int updateLabels(VertexIdx u, Count changeVal);
		VertexIdx getMaxOutNbr(VertexIdx u);

		int insertEdge(eTupleUnWeighted);
		int deleteEdge(eTupleUnWeighted);

		int updateIncPointer(VertexIdx u, VertexIdx v);
		int updateDecPointer(VertexIdx u, VertexIdx v);
		int updateTightInNbrIterator(VertexIdx u, VertexIdx v);
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
		InNbrs.push_back(eachRow);
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

int Graph :: addToPriorityQueue(VertexIdx u, VertexIdx v, Count vVal)
{
	// we need to add/update v in the priority queue of u;
	Count oldVal = nodeToOutdegMap[u][v];
	nodeToOutdegMap[u][v] = vVal;

	outdegToNodeMap[u][oldVal].erase(v);
	outdegToNodeMap[u][vVal].insert(v);

	return 0;
}


int Graph :: removeFromPriorityQueue(VertexIdx u, VertexIdx v)
{
	// remove/decrement v in the priority queue of u;
	Count oldVal = nodeToOutdegMap[u][v];
	Count newVal = oldVal - 1;
	nodeToOutdegMap[u][v] = newVal;

	outdegToNodeMap[u][oldVal].erase(v);
	outdegToNodeMap[u][newVal].insert(v);

	return 0;
}


int Graph :: insertEdge(eTupleUnWeighted e)
{
	// e = u,v -- at this moment an edge isnt directed...
	// hereafter it will be... from here onwards we would orient the edges...
	VertexIdx u = get<0>(e);
	VertexIdx v = get<1>(e);
	VertexIdx w, wPrime;

	if(nodeInDeg[u] >= nodeInDeg[u])
	{
		eTupleUnWeighted eAdd(u,v);
		addDirectedEdge(eAdd);
		w = v;	
	}
	else
	{
		eTupleUnWeighted eAdd(v,u);
		addDirectedEdge(eAdd);
		w = v;
	}

	wPrime = getTightInNbr(w);
	while(wPrime != -1)
	{
		eTupleUnWeighted eFlip(wPrime, w);
		flipDirectedEdge(eFlip);
		w = wPrime;
		wPrime = getTightInNbr(w);
	}

	incrementDu(w);
	return 0;
}


int Graph :: deleteEdge(eTupleUnWeighted e)
{
	// e = u,v directed
	VertexIdx u = get<0>(e);
	VertexIdx v = get<1>(e);
	VertexIdx w, wPrime;

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

	wPrime = getTightOutNbr(w);
	while(wPrime != NullVertexIdx)
	{
		eTupleUnWeighted eFlip(w, wPrime);
		flipDirectedEdge(eFlip);
		w = wPrime;
		wPrime = getTightOutNbr(w);
	}

	decrementDu(w);

	return 0;
}


int Graph :: updateIncPointer(VertexIdx u, VertexIdx v)
{
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

int Graph :: updateDecPointer(VertexIdx u, VertexIdx v)
{
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


int Graph :: updateTightInNbrIterator(VertexIdx u, VertexIdx v)
{
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




int Graph :: addToInNbrs(VertexIdx u, VertexIdx v)
{
	// This is to add v to the in-neighbors of u...
	// add to the list of in-neighbors and update the existence and 
	// address of this new neigbor in the map...
	listOfNeighbors[u].push_back(v);				// adding v to the in-neighbors of u...
	list<VertexIdx>::iterator vsAddress = listOfNeighbors[u].end();
	--vsAddress;									// going to the last element -- which is the new one inserted...
	mapToNeighborsList[u][v] = vsAddress;

	// if this is the first element that is getting added to the list...
	// let the next-pointers point to the first element/neighbor....
	if(listOfNeighbors[u].size() == 1)
	{
		nextPositionIteratorInc[u] = vsAddress;
		nextPositionIteratorDec[u] = vsAddress;
	}

	return 0;
}


int Graph :: removeFromInNbrs(VertexIdx u, VertexIdx v)
{
	// This is to remove v from the in-neighbors of u...

	// check if the current nextIterator is not at the same position as the 
	// node to be removed...
	// if so shift the position of the next-iterator to the next element
	if(mapToNeighborsList[u].find(v) != mapToNeighborsList[u].end())
	{
		if(nextPositionIteratorInc[u] == mapToNeighborsList[u][v])
		{
			updateIncPointer(u, v);
		}

		if(nextPositionIteratorDec[u] == mapToNeighborsList[u][v])
		{
			updateDecPointer(u, v);
		}

		if(nextPositionIteratorTightInNbr[u] == mapToNeighborsList[u][v])
		{
			updateTightInNbrIterator(u, v);
		}

		listOfNeighbors[u].erase(mapToNeighborsList[u][v]);		// removing v from the in-neighbors of u....
		mapToNeighborsList[u].erase(v);							// remove element from the map as well
	}

	return 0;
}


int Graph :: addDirectedEdge(eTupleUnWeighted e)
{
	// e = u,v directed
	VertexIdx u = get<0>(e);
	VertexIdx v = get<1>(e);

	// add u to in-neighbors of v
	InNbrs[v][u] += 1;

	// Note that this is to be done only for the first time when u becomes in-neighbor of v..
	// remember this is a multigraph.... 
	if(InNbrs[v][u] == 1)
	{
		addToInNbrs(v, u);
	}

	// add v to priority queue out-neighbors of u;
	// get the current value of d(v)
	Count dvVal = nodeInDeg[v];

	// update v in the priority queue of u with this new value...
	addToPriorityQueue(u, v, dvVal);

	return 0;
}


int Graph :: removeDirectedEdge(eTupleUnWeighted e)
{
	// e = u,v directed
	VertexIdx u = get<0>(e);
	VertexIdx v = get<1>(e);

	// remove/decrement u from in-neighbors of v
	InNbrs[v][u] -= 1;

	// Note that this is to be done only if the u no more remains an in-neighbor of v..
	// remember this is a multigraph.... 
	if(InNbrs[v][u] == 0)
	{
		removeFromInNbrs(v, u);
	}

	// remove/decrement v from the priority queue out-neighbors of u
	removeFromPriorityQueue(u, v);
	return 0;
}


int Graph :: flipDirectedEdge(eTupleUnWeighted e)
{
	removeDirectedEdge(e);

	VertexIdx u = get<0>(e);
	VertexIdx v = get<1>(e);
	eTupleUnWeighted flippedEdge (v,u);
	addDirectedEdge(e);
	
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

VertexIdx Graph :: getTightInNbr(VertexIdx u)
{
	VertexIdx neighborToReturn = NullVertexIdx;

	list<VertexIdx> :: iterator updateItInc;

	Count start = 0;
	Count maxNumNeighborsToCheck = (Count)(4 * nodeInDeg[u]) / eta;
	
	unordered_map<VertexIdx, int> touchedNeighbors;
	// note that you shouldnt be updating an element multiple times within a same 
	// update... so keep track of the updated neighbors... 
	// every time start with a new/empty list of neighbors that would be updated...
	
	updateItInc = nextPositionIteratorTightInNbr[u];

	while(start < maxNumNeighborsToCheck)
	{
		// access the next neighbor 
		VertexIdx usNextNeighbor = *updateItInc;
		
		if(touchedNeighbors.find(usNextNeighbor) == touchedNeighbors.end())
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
			touchedNeighbors[usNextNeighbor] = 1;

			// check if this current neighbor satisfies the tight in-neighbor condition...
			if(nodeInDeg[usNextNeighbor] <= nodeInDeg[u] - eta/2)
			{
				neighborToReturn = usNextNeighbor;
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

	return neighborToReturn;

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
	string graphFileName = "g1.txt";
	ifstream graphFile;

	string line;
	stringstream ss;
	
	Count n, m;
	VertexIdx e_src, e_dest;
	Weight e_weight, Wmax;
	int insDel;

	double scalingProbParam_c, sparsifyProbParam_c, scalingEps, sparsifyEps;
	double scalingEpsSquared = scalingEps * scalingEps;
	double sparsifyEpsSquared = sparsifyEps * sparsifyEps;

	double scalingProb, minWeightedDensity; 

	cout << "Did I come here...\n";

	graphFile.open(graphFileName, ifstream::in);
	if(graphFile.is_open())
	{
		// if (getline(graphFile, line) != -1)
		if (getline(graphFile, line))
		{
			ss.clear();
			ss.str("");
			ss << line;
			ss >> n >> m >> Wmax;

			// for now, let us say that only the exact number of nodes is only known... 
			// initialize the datastructure required for each node...
			// d(u) 
			// InNbrs -- List of in-nbrs for each node
			// OutNbrs -- Max priority queue for each node u... indexed by d(v)

			Graph G(n); 	// This initializes du, InNbrs, OutNbrs(list of priority queues...)
			cout << "Done with initialization...\n";
			/*

			minWeightedDensity = (Wmax*1.0)/2;

			// scalingProb = (c log n) / (eps^2 * rho_min)  
			scalingProb = (scalingProbParam_c * log2(n)) / (scalingEpsSquared * minWeightedDensity);

			while(getline(graphFile, line))
			{
				ss.clear();
				ss.str("");
				ss << line;
				
				ss >> e_src >> e_dest >> e_weight >> insDel;
				
				// binomial sampling to scale the weight....
				// sparsify -- sampling to decide for each copy of edge to be in the graph or not...
				// maintain sparsified graph -- with edge directions....

				Weight scaledWeight = sampleFromBinomial(e_weight, scalingProb);

				for(Weight w = 0; w < scaledWeight; w++)
				{

				}
			}
			*/
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