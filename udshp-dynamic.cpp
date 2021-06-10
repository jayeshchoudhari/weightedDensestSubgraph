#include "GraphScheduler.hpp"
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <assert.h>
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
	
	double epsUD = std::stod(argv[3]);

    return 0;
}