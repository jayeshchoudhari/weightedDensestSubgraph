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
	Count duplicationFactor = (Count)localAlpha;

	vectorListMap mainEdge2Ids;

	GraphLoader GL(graphFileName);
	Count n = GL.getNumVertices();
	Count numEdgesToBeAdded = GL.getNumOfEdgeAdditionsToBePerformed();
	Count numEdgesToBeRemoved = GL.getNumOfEdgeDeletionsToBePerformed();

	int logn = (int)log2(n);
	int numInstances = logn + 1;
	int maxInstanceId = numInstances - 1;
	EdgeManager EM;

	std::vector<DynamicGraph> DGVecInstance;
	// DynamicGraph DG[numInstances];
	for(int i = 0; i <= maxInstanceId; ++i)
	{
		DGVecInstance.push_back(DynamicGraph(i, n, epsUD, numEdgesToBeAdded, duplicationFactor));
	}

	DynOpManager DOM(maxInstanceId, outFileName);
	DOM.bindGraph(DGVecInstance);
	std::cout<<"Done binding...\n";

	EdgeIdx edgeId = 0;
	EdgeIdx edgeDupId = 0;
	
	auto startTime = std::chrono::system_clock::now();
	Count numOpPerWindow = 0;
	std::cout << "Going for operations.....\n";
	while (GL.has_next())
	{
        EdgeUpdate edge_up = GL.next_update();
        // int edge_id = 0;
        if (!edge_up.is_report) 
		{
			numOpPerWindow += 1;
			///////////// INSERTION /////////////////////////
			// std::vector<VertexIdx> eVec = edge_up.vertices;
            if (edge_up.is_add) 
			{
				// std::cout << "insertion -- " << edgeId << "\n";
				mainEdge2Ids[edge_up.vertices].push_back(edgeId);
				std::vector<EdgeIdx> edgeDuplicatorIds = EM.getEdgeIdsAfterDuplication(edge_up.vertices, duplicationFactor);
                DOM.addEdge(edge_up.vertices, edgeDuplicatorIds, EM);
				edgeId += 1;
            }
			else
			{
				///////////////// DELETION /////////////////////////
				if(mainEdge2Ids.find(edge_up.vertices) != mainEdge2Ids.end())
				{
					// std::cout << "Deletion --- " << std::endl;
					std::vector<EdgeIdx> edgeDuplicatorIds = EM.retrieveDuplicatedEdgeIds(edge_up.vertices, duplicationFactor);
					// std::cout << "Got the edge Ids --- " << std::endl;
					// for(int i = 0; i < edge_up.vertices.size(); i++)
					// 	std::cout << edge_up.vertices[i] << " ";
					// std::cout << "\n";
					DOM.removeEdge(edge_up.vertices, edgeDuplicatorIds, EM);
					// std::cout << "Edge Deleted --- " << std::endl;
					// EM.removeEdgeFromMemory(edgeDuplicatorIds);
					// std::cout << "Edges Removed from memory --- " << std::endl;
					// remove edge from the main set of edgeIds 
					mainEdge2Ids[edge_up.vertices].pop_front();
					if(mainEdge2Ids[edge_up.vertices].size() == 0)
					{
						mainEdge2Ids.erase(edge_up.vertices);
					}
				}
				else
				{
					std::cout << "This edge was never added before... \n";
				}
            }
        }
		else
		{
			std::cout << "Reporting...\n";
			DOM.addPendingEdgesToActiveInstance(EM);
			auto endTime = std::chrono::system_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
			double micros = elapsed.count();
			DOM.getDensityEstimate(EM, localAlpha, mainEdge2Ids, micros, numOpPerWindow, edge_up.report_label);
			auto lastTime = std::chrono::system_clock::now();
			std::cout << "Done with first report.. TIME -- " << micros/numOpPerWindow << "\n"; 
			numOpPerWindow = 0;
			startTime = lastTime;
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