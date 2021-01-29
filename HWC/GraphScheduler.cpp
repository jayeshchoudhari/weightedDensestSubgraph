/*
 *  GraphScheduler.cpp
 *  Modified by sghu@cs.hku.hk
 *  Created on: Oct 22, 2014
 *      Author: aepasto
 */

#include "GraphScheduler.hpp"
#include "Utility.hpp"

#include <string>
#include <cassert>
#include <iostream>
#include <functional>
using namespace std;

GraphScheduler::~GraphScheduler() {
	this->file_stream_.close();
}

GraphScheduler::GraphScheduler(const string& file_name) {

	//file_stream_.open(file_name.c_str(), ios_base::in);
    file_stream_.open(file_name.c_str(), ios::in);
	if (file_stream_.is_open())
	    std::cout<< "open file" <<std::endl;
	retrieve_all_edges();
	position_queue_ = 0;
	add_count = 0;
	remove_count = 0;
}

void GraphScheduler::retrieve_all_edges() {

	std::string delimiter = " ";
	string line;
	vector<string> tokens;


	while (getline(file_stream_, line)) {
		tokens.clear();

		size_t pos = 0;
		std::string token;
		while (pos != line.npos) {
			pos = line.find(delimiter);
			token = line.substr(0, pos);
			tokens.push_back(token);
			line.erase(0, pos + delimiter.length());
		}

		//assert(tokens.size() >= 2);


        EdgeUpdate next_edge;

		if (tokens[0][0] == '+') {
			next_edge.is_add = true;
            next_edge.is_report = false;
		} else if (tokens[0][0] == '-') {
			next_edge.is_add = false;
            next_edge.is_report = false;
		}
		else if (tokens[0] == "report") {
            next_edge.is_add = false;
		    next_edge.is_report = true;
		    next_edge.report_label = tokens[1];
        }
		else {
//		    std::cout << "Unknown token, ignoring the line" << std::endl;
            continue;
//		    assert(false);
		}
        if (!next_edge.is_report) {
            for (size_t i = 1; i < tokens.size(); i++)
                next_edge.vertices.push_back(atoi(tokens[i].c_str()));

            if (next_edge.is_add) {
                next_edge.timestamp = next_edge.vertices.back();
                next_edge.vertices.pop_back();
            }


            sort(next_edge.vertices.begin(), next_edge.vertices.end());

            next_edge.vertices.shrink_to_fit();
        }

        edge_queue_.push_back(next_edge);
		

        const int MOD = 1000000; // 100000000
		if (edge_queue_.size() % MOD == MOD-1) {
			cerr << "READ #: " << edge_queue_.size() + 1 << endl;
		}
	}

	// Only in C++11
	edge_queue_.shrink_to_fit();
	//edge_queue_no_time_.shrink_to_fit();
}

EdgeUpdate GraphScheduler::next_update() {

	EdgeUpdate edge_queue;
	edge_queue = edge_queue_[position_queue_];

	++position_queue_;

	if (!edge_queue.is_report) {
        if (edge_queue.is_add) {
            ++add_count;
        } else {
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
