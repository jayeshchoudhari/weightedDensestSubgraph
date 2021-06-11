#include "./include/namespace.h"
#include "./include/GraphLoader.h"
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

	GraphLoader GL(graphFileName);

	Count n = GL.getNumVertices();


	int logn = (int)log2(n);
	int numInstances = logn + 1;
	int maxInstanceId = numInstances - 1;


	DynamicGraph Ti[numInstances];
	for(int i = 0; i <= maxInstanceId; ++i)
	{
		Ti[i].initializeVariables(i, n, epsUD);
	}

    return 0;
}