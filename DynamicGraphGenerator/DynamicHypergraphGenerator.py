from typing import NamedTuple


class DynInputParam(NamedTuple):
    duplicate: bool
    incremental: bool
    out_f_path: str
    time_window: int = 1
    sliding_window_factor: int = 1  # Will only be used if the timestamp semantics is unclear


# Read from a hypergraph file and store it into a list of edges
def read_edge_file(f_path):
    edge_list = []
    with open(f_path, 'rt') as f:
        for line in f:
            if line[0] != '%':
                edge = line.split()
                time = int(edge[-1])
                edge = [int(item) for item in edge[:-1]]
                edge_list.append([edge, time])
    return edge_list


# A method to verify that there are no duplicate hyperedges
# at any point in the edgeList
def check_duplicate(edge_list):
    active_edge_set = {}
    for edge in edge_list:
        if edge[0] == '+':
            n_edge = edge[1:-1]
            n_edge = [int(vertex) for vertex in n_edge]
            n_edge.sort()
            if str(n_edge) in active_edge_set:
                print("Duplicate edge:" + str(edge))
                return False
            else:
                active_edge_set[str(n_edge)] = 1
        elif edge[0] == '-':
            n_edge = edge[1:]
            n_edge = [int(vertex) for vertex in n_edge]
            n_edge.sort()
            if str(n_edge) not in active_edge_set:
                print("Removing non existent edge:" + str(edge))
                return False
            else:
                del active_edge_set[str(n_edge)]
    return True


# Using the hyperedges construct a dynamic hypergraph
# In order construct the dynamic hypergraph, we consider a
# sliding window model. The timestamp for each edge is used
# to decide the sliding window. It is assumed that the time_window and
# the timestamps corresponding to each hyperedge is scaled appropriately.
# It is suggestive to think of timestamp as year or quarter+year
# or month+year (converted  to an integer)
def make_dynamic_hg_by_timestamp(hg_edge_list: list, dyn_input_param: DynInputParam):

    # Each element edge_list is of type [edge], timestamp
    # Each entry in the edge is an integer
    # hg_edge_list is sorted by timestamp
    hg_edge_list.sort(key=lambda k:k[1])

    # Hash the authors to integer values in the sorted list
    index = 0
    dyn_edge_list = []
    author_map = {}
    active_edge_set = {}
    edges_by_timestamp = {}
    earliest_timestamp = hg_edge_list[0][1]

    # The following two maps are for book-keeping purpose
    # They maintain the number of hyperedges added
    # and removed at each timestamp
    edge_added_by_year = {}
    edge_removed_by_year = {}

    max_card = 0
    max_edge_count = 0

    for vertices, timestamp in hg_edge_list:

        if timestamp not in edge_added_by_year:
            edge_added_by_year[timestamp] = 0
        dup = False

        if timestamp - earliest_timestamp >= dyn_input_param.time_window:

            # Add a report statement before every new time-stamp
            # The current graph is a sliding of size dyn_input_param.time_window
            # ending just before the current time_stamp
            dyn_edge_list.append(["report "+str(timestamp-1)])
            max_edge_count = max(max_edge_count, len(active_edge_set))

            if timestamp not in edge_removed_by_year:
                edge_removed_by_year[timestamp] = 0

            # For fully-dynamic input graph, remove hyperedges
            # The hyperedges to be removed are exactly the ones
            # with earliest_timestamp
            # We remove all the edges first and then subsequently
            # add edges with timestamp as the current timestamps
            if not dyn_input_param.incremental:
                for edge in edges_by_timestamp[earliest_timestamp]:
                    if not dyn_input_param.duplicate:
                        active_edge_set.pop(str(edge))
                    # add the deletion symbol at the beginning
                    mod_edge = ['-']+edge
                    dyn_edge_list.append(mod_edge)
                    edge_removed_by_year[timestamp] += 1

            earliest_timestamp += 1  # We are assuming presence of consecutive timestaamps

        # We add the current hyperedge to the graph
        edge = []
        for vertex in vertices:
            if vertex not in author_map:
                author_map[vertex] = index
                index += 1
            edge.append(author_map[vertex])
            edge.sort()  # sorting this hyperedge ensure no duplicate edges, if required

        if str(edge) not in active_edge_set:
            active_edge_set[str(edge)] = 1
        else:
            dup = True

        if dyn_input_param.duplicate or not dup:
            if timestamp not in edges_by_timestamp:
                edges_by_timestamp[timestamp] = [edge]
            else:
                edges_by_timestamp[timestamp].append(edge)
            max_card = max(max_card, len(edge))

            # Insert the edge along with the timestamp
            mod_edge = ['+']+edge+[timestamp]
            dyn_edge_list.append(mod_edge)

            edge_added_by_year[timestamp] += 1

    # Finally add one more report statement to the end of the graph stream
    dyn_edge_list.append(["report " + str(hg_edge_list[-1][1])])

    # no of distinct vertices in the entire hypergraph
    n = index
    # If the hypergraph is not supposed to have any duplcaites,
    # then ensure that the edge_list is duplicate free in a sliding window
    if dyn_input_param.duplicate or check_duplicate(dyn_edge_list):
        print("Passed sanity check")

    # Now write the hypergraph to a txt file
    write_hg_file(dyn_edge_list, n, max_card, max_edge_count, dyn_input_param)

    return dyn_edge_list, n, max_card, max_edge_count


# This function will be used if the semantics of the timestamp is not clear
# In that case, we will simply divide the total number of hyperedges
# by sliding_window_factor. Then we will add part 1, then remove part 1
# and part 2 and so on. Here, the addition and the removal of the edges
# are interleaved. For intance, let's say sliding_window_factor posits that
# each part is of size k. Then, we first add first k hyperedges to the graph,
# then remove 1st hyperedge and add k+1 th hyperedge and so on.
def make_dynamic_hg_by_size(hg_edge_list:list, dyn_input_param:DynInputParam):

    # Each element edge_list is of type [edge], timestamp
    # Each entry in the edge is an integer
    # hg_edge_list is sorted by timestamp
    hg_edge_list.sort(key=lambda k:k[1])

    max_card = 0
    max_edge_count = 0
    dyn_edge_list = []
    vertex_set = set()
    active_edge_set = {}

    start_time = hg_edge_list[0][1]
    end_time = hg_edge_list[-1][1]

    sliding_window_size = 0
    report_freq = 0
    sliding_window_size = int(len(hg_edge_list) / dyn_input_param.sliding_window_factor)
    report_freq = sliding_window_size
    print (hg_edge_list[-1][1], hg_edge_list[0][1], len(hg_edge_list), report_freq)

    left_window = 0
    for i,item in enumerate(hg_edge_list):
        edge,time = item
        dup = False
        if i % report_freq == report_freq-1:
            dyn_edge_list.append(["report "+str(time-1)])
            print ("report:"+str(i)+'time:'+str(time))
            max_edge_count = max(max_edge_count,len(active_edge_set))

        if not dyn_input_param.incremental and i > sliding_window_size:
            delete_edge =hg_edge_list[left_window][0]
            delete_edge.sort()
            if not dyn_input_param.duplicate:
                if str(delete_edge) in active_edge_set:
                    active_edge_set.pop(str(delete_edge))
                    # add the deletion symbol at the beginning
                    mod_edge = ['-']+delete_edge
                    dyn_edge_list.append(mod_edge)
            left_window += 1

        edge.sort()
        vertex_set.update(edge)
        if str(edge) not in active_edge_set:
            active_edge_set[str(edge)] = 1
        else:
            dup = True

        if dyn_input_param.duplicate or not dup:
            max_card = max(max_card,len(edge))
            mod_edge = ['+']+edge+[time]
            dyn_edge_list.append(mod_edge)

    dyn_edge_list.append(["report " + str(hg_edge_list[-1][1])])
    # Here n does not repsent the length of the true vertex set and
    # should be used as such
    n = len(vertex_set)
    print(max_edge_count,report_freq)

    # If the hypergraph is not supposed to have any duplcaites,
    # then ensure that the edge_list is duplicate free in a sliding window
    if dyn_input_param.duplicate or check_duplicate(dyn_edge_list):
        print("Passed sanity check")

    # Now write the hypergraph to a txt file
    write_hg_file(dyn_edge_list, n, max_card, max_edge_count, dyn_input_param)

    return dyn_edge_list,n,max_card,max_edge_count


# Write the dynamic hypergraph to a file alongwith other
# header information
def write_hg_file(edge_list, n, max_card, max_edge_count, dyn_input_param):

    with open(dyn_input_param.out_f_path, 'wt') as outFile:
        # Write header information
        outFile.write('% ' + str(dyn_input_param) + '\n')
        outFile.write('% n = ' + str(n) +
                      ', max cardinality = ' + str(max_card) +
                      ',max_edge_count = ' + str(max_edge_count) + '\n')
        for edge in edge_list:
            # Remove single author paper here
            str_edge = [str(element) for element in edge]
            outFile.write(' '.join(str_edge))
            outFile.write('\n')


if __name__ == '__main__':

    file_id = 'tags-stack-overflow'

    dir_path = 'data/'+file_id+'/'
    # file_id = file_id + '.all.1985'  # Only for dblp dataset
    input_file_path = dir_path + file_id + '-hg.txt'

    edge_list = read_edge_file(input_file_path)

    # time_window dictates the no of years in consideration for the hypergraph
    # The desnity will be reported for each sliding window of time_window many duration
    # For eg, consider time_window = 5 and incremental = false
    # Assume we are looking at the dblp dataset with timestamp represeting a year
    # first the papers from 1985-1989 will be added. There will a report command
    # Then, papers from 1985 will be removed and 1990 will be added
    time_windows = 1
    window_factor = 1 #this to be used only if the data is to be divided equally

    # Does the hypergraph support having multiple hyperedges?
    duplicate = False

    # If incremental is set to yes, then only hyperedge insertion will occur and
    # no hyperedge deletion will be done
    incremental = True

    dyn = 'D'
    if incremental:
        dyn = 'I'

    hg_f_name = dir_path + file_id + '.hg.'+dyn+'.'+str(time_windows)+'.txt'
    input_param = DynInputParam(duplicate=duplicate,
                                incremental=incremental,
                                out_f_path=hg_f_name,
                                time_window=time_windows,
                                sliding_window_factor=window_factor)  # This parameter will not be used
    make_dynamic_hg_by_timestamp(edge_list, input_param)

    # hg_f_name = dir_path + file_id + '.hg.' + dyn + '.' + str(window_factor) + '.txt'
    # input_param = DynInputParam(duplicate=duplicate,
    #                             incremental=incremental,
    #                             out_f_path=hg_f_name,
    #                             time_window=time_windows, # This parameter will not be used
    #                             sliding_window_factor=window_factor)
    # make_dynamic_hg_by_size(edge_list, input_param)
