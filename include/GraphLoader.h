#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

#include "namespace.h"

// using namespace std;

class GraphLoader
{
    private:
        Count add_count;
        Count remove_count;
        Count loadAddCount, loadRemoveCount;
        unsigned int position_queue_;
        std::ifstream file_stream_;
        std::vector<EdgeUpdate> edge_queue_;
        Count numVertices;
        unsigned int maxRank;

    public:
        GraphLoader(const std::string& file_name);
        // virtual ~GraphScheduler();

        int loadAllEdges();

        EdgeUpdate next_update();

        Count getNumVertices()
        {
            return numVertices;
        }
        
        Count getNumOfEdgeAdditionsToBePerformed()
        {
            return loadAddCount;
        }

        Count getNumOfEdgeDeletionsToBePerformed()
        {
            return loadRemoveCount;
        }

        int getMaxEdgeRank()
        {
            return maxRank;
        }

        inline bool has_next()
        {
            return (unsigned int) position_queue_ < (unsigned int) edge_queue_.size();
        }
};



#endif  // GRAPH_LOADER_H