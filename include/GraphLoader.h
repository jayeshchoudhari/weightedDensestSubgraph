#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

typedef struct EdgeUpdate
{
    std::vector<int> vertices;
    int timestamp;
	bool is_add;
	bool is_report;
	std::string report_label;
} EdgeUpdate;


class GraphLoader
{
    private:
        int add_count;
        int remove_count;
        int loadAllEdges();
        ifstream file_stream_;
        unsigned int position_queue_;
        vector<EdgeUpdate> edge_queue_;

    public:
        GraphScheduler(const string& file_name);
        // virtual ~GraphScheduler();

        EdgeUpdate next_update();
        inline bool has_next()
        {
            return (unsigned int) position_queue_ < (unsigned int) edge_queue_.size();
        }
};



#endif  // GRAPH_LOADER_H