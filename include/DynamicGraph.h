#ifndef DYNAMIC_GRAPH_H
#define DYNAMIC_GRAPH_H

#include "namespace.h"
#include "EdgeManager.h"
#include "GraphLoader.h"
#include <iostream>

class DynamicGraph 
{
    private:
        Count nVertices; 								// number of vertices in the graph
        Count nEdges;     								// number of edges in this list
        double alpha;									// multiplicity
        double epsVal;									// error factor
		double eta;										// relaxation
 		double rhoEst;									// estimated density


        static std::multimap<std::vector<VertexIdx>, EdgeIdx> edgeMap;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        static std::vector <std::vector<VertexIdx>> edgeList;

        // above could be static members shared among all the instances of the Graph...
        // each instance for a guess of density...
	
        // std::unordered_map<EdgeIdx, int> pendingListOfEdges;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
        std::unordered_map<VertexIdx, std::set<EdgeIdx>> pendingForNode;
        std::unordered_map<EdgeIdx, std::set<VertexIdx>> nodesWithPendingEdge;

 		

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
		std::unordered_map<EdgeIdx, std::vector<VertexIdx>> pendingListOfEdges;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
		// std::unordered_map<EdgeIdx, int> pendingListOfEdges;	// Map of edges -- to check mostly if the edge already exists or not...and to keep an Id for each edge...
		DynamicGraph(int i, Count nv, float epsUD);
		
    	// static int initializeStaticMembers(Count nv, double decEps);
    	int initializeVariables(int i, Count nv, float epsVal);

		int insertEdge(edgeVector &e, EdgeIdx eId, EdgeManager &EM);
		// int deleteEdge(edgeVector &e, EdgeIdx eId, EdgeManager &EM);
		VertexIdx deleteEdge(edgeVector &e, EdgeIdx eId, EdgeManager &EM);

    	int showInstanceVariables();

    	double getRhoEst();

    	int addEdgeToPendingList(edgeVector &e, EdgeIdx eId);
    	int checkEdgeExistenceInPendingList(EdgeIdx eId);
		int removeEdgeFromPendingList(edgeVector &e, EdgeIdx eId);
		std::pair<EdgeIdx, edgeVector> getPendingEdgeForLastVertex(VertexIdx lv);
		// EdgeIdx getPendingEdgeForLastVertex(VertexIdx lv);

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

		int addDirectedEdgeToInOutNbrs(edgeVector &currentEdge, EdgeIdx eId, VertexIdx newHeadNode);
		int removeDirectedEdgeFromInOutNbrs(edgeVector &currentEdge, EdgeIdx eId, VertexIdx headNode);
		int flipDirectedEdge(EdgeIdx eId, VertexIdx oldHeadNode, VertexIdx newHeadNode, EdgeManager &EM);
		
		int addToPriorityQueue(VertexIdx, VertexIdx, Count, EdgeIdx);
		int removeFromPriorityQueue(VertexIdx, VertexIdx, EdgeIdx);
		
		int updateNextNeighbors(VertexIdx u, Count newDuVal, int incOrDec, EdgeManager &EM);

		int incrementDu(VertexIdx u, EdgeManager &EM);
		int decrementDu(VertexIdx u, EdgeManager &EM);

		std::pair<VertexIdx, EdgeIdx> getTightInNbr(VertexIdx v, EdgeManager &EM);
		VertexIdx getMinDegreeVertexInE(edgeVector &currentEdge);
		Count getMinLoadInE(EdgeIdx eId, EdgeManager &EM);

		EdgeIdx getTightOutNbr(VertexIdx);
		EdgeIdx getMaxOutNbr(VertexIdx u);

		int updateLabels(VertexIdx u, Count changeVal);

		

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

		double getDensityOfInducedSubgraph(std::set<VertexIdx> denseSubgraphNodes, vectorListMap &mainEdge2Ids, Count duplicationFactor);
		std::pair<double, unsigned int> getMaxPartitionDensity(vectorListMap &mainEdge2Ids, Count duplicationFactor);

		unsigned int getPendingCount();
		int insertListOfPendingEdges(EdgeManager &EM);

		int showPQs();
};

#endif  // GRAPH_H