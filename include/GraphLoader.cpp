#include "namespace.h"
#include "GraphLoader.h"
#include <iostream>
#include <string>
#include <cassert>
#include <stdlib.h>

// using namespace std;

GraphLoader::GraphLoader(const std::string& file_name) 
{
	//file_stream_.open(file_name.c_str(), ios_base::in);
    file_stream_.open(file_name.c_str(), std::ios::in);
	if (file_stream_.is_open())
    {
	    std::cout<< "open file" <<std::endl;
        loadAddCount = 0;
        loadRemoveCount = 0;
        loadAllEdges();
        position_queue_ = 0;
        add_count = 0;
        remove_count = 0;
    }
    else
    {
        std::cout << "Cannot open input file...\n";
		exit(0);
    }
}


int GraphLoader::loadAllEdges() 
{
	std::string delimiter = " ";
	std::string line;
	std::vector<std::string> tokens;
    
    int lineno = 0;

    while (getline(file_stream_, line)) 
    {
        size_t pos = 0;
        char *eptr;
        if(lineno == 0)
        {
            std::string token_n, token_mk;
            pos = line.find(delimiter);
            token_n = line.substr(0, pos);
            token_mk = line.substr(pos, line.length());
            numVertices = atoi(token_n.c_str());
            maxRank = atoi(token_mk.c_str());
        }
        else
        {
            tokens.clear();
            std::string token;
            while (pos != line.npos) 
            {
                pos = line.find(delimiter);
                token = line.substr(0, pos);
                tokens.push_back(token);
                line.erase(0, pos + delimiter.length());
            }

            //assert(tokens.size() >= 2);
            EdgeUpdate next_edge;

            if (tokens[0][0] == '+') 
            {
                next_edge.is_add = true;
                next_edge.is_report = false;
                loadAddCount += 1;
            } 
            else if (tokens[0][0] == '-') 
            {
                next_edge.is_add = false;
                next_edge.is_report = false;
                loadRemoveCount -= 1;
            }
            else if (tokens[0] == "=")
            {
                next_edge.is_add = false;
                next_edge.is_report = true;
                next_edge.report_label = tokens[1];
            }
            else 
            {
    //		    std::cout << "Unknown token, ignoring the line" << std::endl;
                continue;
    //		    assert(false);
            }
            if (!next_edge.is_report) 
            {
                for (size_t i = 1; i < tokens.size(); i++)
                {
                    // next_edge.vertices.push_back(atoi(tokens[i].c_str()));
                    next_edge.vertices.push_back(strtol(tokens[i].c_str(), &eptr, 10));
                }
                if (next_edge.is_add) 
                {
                    next_edge.timestamp = next_edge.vertices.back();
                    next_edge.vertices.pop_back();
                }

                sort(next_edge.vertices.begin(), next_edge.vertices.end());

                next_edge.vertices.shrink_to_fit();
            }

            edge_queue_.push_back(next_edge);
            
            const int MOD = 1000000; // 100000000
            if (edge_queue_.size() % MOD == MOD-1) 
            {
                std::cerr << "READ #: " << edge_queue_.size() + 1 << std::endl;
            }
        }
    
        lineno += 1;
    }

    // Only in C++11
    edge_queue_.shrink_to_fit();
    //edge_queue_no_time_.shrink_to_fit();


    file_stream_.close();

    return 0;
}


EdgeUpdate GraphLoader::next_update() 
{
	EdgeUpdate edge_queue;
	edge_queue = edge_queue_[position_queue_];

	++position_queue_;

	if (!edge_queue.is_report) 
    {
        if (edge_queue.is_add) 
        {
            ++add_count;
        } 
        else 
        {
            ++remove_count;
        }
	}
	/*if (edge_queue.is_add && add_count % 100000 == 99999) {
		cerr << "ADD #: " << add_count + 1 << endl;
		cerr.flush();
	} else if (!edge_queue.is_add && remove_count % 100000 == 99999) {
		cerr << "REM #: " << remove_count + 1 << endl;
		cerr.flush();
	}*/
	return edge_queue;
}