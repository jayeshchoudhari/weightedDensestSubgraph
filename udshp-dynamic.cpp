#include "./include/namespace.h"
#include "./include/GraphLoader.h"
#include "./include/EdgeManager.h"
#include "./include/DynOpManager.h"
#include "./include/DynamicGraph.h"
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>


int main(int argc, char** argv)
{
	if (argc != 4) 
	{
        std::cerr << "ERROR -- Requires 3 parameters--: 1) Input Graph filename; 2) Output filename; and 3) Epsilon (error parameter);\n";
        exit(1);
    }

	std::string graphFileName = argv[1];

	std::ofstream outFile;
	std::string outFileName = argv[2];
	
	float epsUD = std::stod(argv[3]);
	double localAlpha = round(1.0/pow(epsUD,2));
	int duplicationFactor = (int)localAlpha;
	// vector<vector<VertexIdx>> mainEdgePool;

	std::unordered_map<std::vector<VertexIdx>, std::list<EdgeIdx>, vectorHash> mainEdge2Ids;

	GraphLoader GL(graphFileName);
	Count n = GL.getNumVertices();

	int logn = (int)log2(n);
	int numInstances = logn + 1;
	int maxInstanceId = numInstances - 1;
	EdgeManager EM;

	std::vector<DynamicGraph> DGVecInstance;
	// DynamicGraph DG[numInstances];
	for(int i = 0; i <= maxInstanceId; ++i)
	{
		DGVecInstance.push_back(DynamicGraph(i, n, epsUD));
	}

	DynOpManager DOM(maxInstanceId);
	DOM.bindGraph(DGVecInstance);

	EdgeIdx edgeId = 0;
	EdgeIdx edgeDupId = 0;

	while (GL.has_next())
	{
        EdgeUpdate edge_up = GL.next_update();
        // int edge_id = 0;
        if (!edge_up.is_report) 
		{
			// std::vector<VertexIdx> eVec = edge_up.vertices;
            if (edge_up.is_add) 
			{
				mainEdge2Ids[edge_up.vertices].push_back(edgeId);
				std::vector<EdgeIdx> edgeDuplicatorIds = EM.getEdgeIdsAfterDuplication(edge_up.vertices, duplicationFactor);
                DOM.addEdge(edge_up.vertices, edgeDuplicatorIds, EM);
				edgeId += 1;
            }
			else
			{
				if(mainEdge2Ids.find(edge_up.vertices) != mainEdge2Ids.end())
				{
					std::vector<EdgeIdx> edgeDuplicatorIds = EM.retrieveDuplicatedEdgeIds(edge_up.vertices, duplicationFactor);
					// edge_id = DOM.removeEdge(edge_up.vertices);
				}
            }
        }
        // assert(edge_id < std::numeric_limits<int>::max() );

        // double upper_bound = h.maxEdgeCardinality() * ads.beta() * (1.0 + ads.epsilon());

        // if (edge_id != -1)
		// {
        //     stats.exec_op(edge_up.is_add, edge_up.is_report, ads.subgraphSize(),
        //             ads.density(), upper_bound, edge_up.timestamp, edge_up.report_label);
        // }
    }

    return 0;
}