#include "include/Graph.h"
#include <bits/stdc++.h>

using namespace std;
// using namespace rwNamespace;

class Graph 
{
    private:
        Count nVertices; 								//number of vertices in the graph
        Count nEdges;     								//number of edges in this list
        vector <VertexIdx> srcs;      					//array of source vertices
        vector <VertexIdx> dsts;      					//array of destination vertices
        vector <VertexIdx> nodeList;
        vector <vector<VertexIdx>> adjList;				// adj List
        vector <vector<Weight>> adjMatrix;				//adj Matrix
        vector <vector<Weight>> scaledAdjMatrix;		//scaled adj Matrix
        vector <vector<Weight>> sparsedAdjMatrix;		//sparsed adj Matrix
        vector <vector<VertexIdx>> edgeList;
        unordered_map <VertexIdx, Count> nodeDeg;

        /* Count Lists... */
        // vector<int> g3s(3);
        // vector <Count> g3s = vector<Count>(3);

    public:

    	Graph(vector <VertexIdx> &v, vector<ePair> &edList);		//takes in list/vector of nodes and vector of pairs which are edges (undirected)

		Graph(int nv, int ne);		//gets number of vertices and edges... edges might change -- but this is just about the file... 

    	int printGraphDetails();

  		// int CountG3s();
		// ePair lStepRandomWalk(Count, VertexIdx);

		Count getNumVertices();
		Count getNumEdges();
        // //Make a copy
        // Graph copy() const;

        // //For debugging.  Not for serialization.
        // void print(FILE *f = stdout) const;
};


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
	for(nIt = nodeDeg.begin(); nIt != nodeDeg.end(); nIt++)
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


int Graph :: checkEdgeExistence(e)
{
	// check in the adjacency list... if the pair exists....
	return flag;
}

// // print adjacency list representation of graph
// void printGraph(Graph const& graph, int N)
// {
// 	for (int i = 0; i < N; i++)
// 	{
// 		// print current vertex number
// 		cout << i << " --> ";

// 		// print all neighboring vertices of vertex i
// 		for (int v : graph.adjList[i])
// 			cout << v << " ";
// 		cout << endl;
// 	}
// }


Graph :: Graph(vector <VertexIdx> &v, vector<ePair> &e)
// Graph::Graph(vector <VertexIdx> &v, vector<ePair> &edList)
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

	// vector<ePair>::iterator edgeIt;

	for(const ePair &edgeIt : e)
	{
		VertexIdx src = edgeIt.first;
		VertexIdx dest = edgeIt.second;

		// populate adjacency list
		adjList[src].push_back(dest);
		adjList[dest].push_back(src);

		// populate list of edges
		vector <VertexIdx> tempEdge;
		tempEdge.push_back(src);
		tempEdge.push_back(dest);
		edgeList.push_back(tempEdge);

		// populate node degress
		nodeDeg[src] += 1;
		nodeDeg[dest] += 1;
	}
}

Graph :: Graph(int nv, int ne)
{
	nVertices = nv;

	// initialize adj matrix to 0...
	// each entry would be positive > 0 if an edge has some positive weight...
	for(Count i = 0; i < nv, i++)
	{
		vector <Weight> eachRow;
		for(Count j = 0; j < nv, j++)
		{
			eachRow.push_back(0);
		}
		adjMatrix.push_back(eachRow);
		scaledAdjMatrix.push_back(eachRow);
		sparsedAdjMatrix.push_back(eachRow);
	}
}


int Graph :: CountG3s()
{
	vector <VertexIdx> :: iterator neighborIt;
	vector <VertexIdx> neighbors;

	unordered_map<long long, Count> tri_e;

	for(int e = 0; e < nEdges; e++)
	{
		Count star_u = 0, star_v = 0;

		unordered_map<VertexIdx, bool> X;

		VertexIdx u = edgeList[e][0];
		VertexIdx v = edgeList[e][1];

		// Checking neighbors of u...
		neighbors = adjList[u];
		for(neighborIt = neighbors.begin(); neighborIt != neighbors.end(); neighborIt++)
		{
			VertexIdx w = *neighborIt;
			if (w == v) continue;
			star_u += 1;
			X[w] = 1;
		}

		// Checking neighbors of v now..
		neighbors = adjList[v];
		for(neighborIt = neighbors.begin(); neighborIt != neighbors.end(); neighborIt++)
		{
			VertexIdx w = *neighborIt;
			if (w == u) continue;

			if (X[w] == 1)
			{
				tri_e[e] += 1;
				star_u -= 1;
			}
			else
			{
				star_v += 1;
			}
		}

		g3s[0] += tri_e[e];
		g3s[1] += (star_u + star_v);
		X.clear();
	}

	cout << g3s[0] << " " << g3s[1] << "\n"; 


	return 0;
}


int sampleFromBinomial(int wt, double p)
{
	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// default_random_engine generator (seed);
  	random_device rd;
    mt19937 gen(rd());
    // perform 4 trials, each succeeds 1 in 2 times
    // std::binomial_distribution<> d(4, 0.5);

  	std::binomial_distribution<int> distribution (n, p);
  	return distribution(gen);
}


int processEdge()
{

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
	ifstream graphFile;

	string line;
	stringstream ss;
	
	Count n, m;
	VertexIdx e_src, e_dest;
	Weight e_weight; Wmax;
	int insDel;

	double scalingProbParam_c, sparsifyProbParam_c, scalingEps, sparsifyEps;
	double scalingEpsSquared = scalingEps * scalingEps;
	double sparsifyEpsSquared = sparsifyEps * sparsifyEps;

	double scalingProb, minWeightedDensity; 

	graphFile.open(graphFileName);

	if(graphFile.is_open())
	{
		if (getline(graphFile, line) != -1)
		{
			ss.clear();
			ss.str("");
			ss << line;
			ss >> n >> m >> Wmax;

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
		}
		else
		{
			cout << "First line is empty line in the file...\n";
			exit(0);
		}
		
	}



	Graph G(v, e);

	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	G.printGraphDetails();

	G.CountG3s();


	return 0;
}