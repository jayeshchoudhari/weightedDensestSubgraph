/*
 *  GraphScheduler.h
 *  Modified by sghu@cs.hku.hk
 *  Created on: Oct 22, 2014
 *      Author: aepasto
 */

#ifndef GRAPHSCHEDULER_H_
#define GRAPHSCHEDULER_H_

#include <algorithm>
#include <vector>
#include <queue>
#include <fstream>
#include <string>
using namespace std;

enum Update {
	ADD, REM
};

typedef struct EdgeUpdate {
	std::vector<int> vertices;
    int timestamp;
	bool is_add;
	bool is_report;
	std::string report_label;
} EdgeUpdate;


class GraphScheduler {
public:
	GraphScheduler(const string& file_name);
	virtual ~GraphScheduler();

	EdgeUpdate next_update();
	inline bool has_next() {
	   return (unsigned int) position_queue_
			  < (unsigned int) edge_queue_.size();
	}

private:
	int add_count;
	int remove_count;
	void retrieve_all_edges();
	ifstream file_stream_;
	vector<EdgeUpdate> edge_queue_;
	unsigned int position_queue_;
};

#endif /* GRAPHSCHEDULER_H_ */
