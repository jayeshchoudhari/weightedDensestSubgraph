from typing import NamedTuple
from DynamicGraphGenerator import DBLPjson2Hypergraph as dbLPH


class DynInputParam(NamedTuple):
    duplicate: bool
    incremental: bool
    sliding_window_factor: int
    out_f_path: str

def make_dynamic_hg(edge_list, dyn_input_param):

    # Each element edge_list is of type [edge], time_label
    # Each entry in the edge is an integer
    # paper_list is sorted by time_label
    edge_list.sort(key=lambda k:k[1])

    max_card = 0
    max_edge_count = 0
    edgeList = []
    vertex_set = set()
    active_edge_set = {}

    start_time = edge_list[0][1]
    end_time = edge_list[-1][1]

    sliding_window_size = 0
    report_freq = 0
    if dyn_input_param.sliding_window_factor == -1:
        # sliding_window_size = int((end_time-start_time)/200)
        # report_freq = int((end_time-start_time)/200)
        sliding_window_size = int(len(edge_list)/20)
        report_freq = int(len(edge_list)/20)
        print (edge_list[-1][1], edge_list[0][1], len(edge_list),report_freq)
    else:
        sliding_window_size = int(len(edge_list)/dyn_input_param.sliding_window_factor)
        report_freq = int(len(edge_list)/dyn_input_param.sliding_window_factor)
        print (edge_list[-1][1], edge_list[0][1], len(edge_list),report_freq)

    left_window = 0
    for i,item in enumerate(edge_list):
        edge,time = item
        dup = False
        if i % report_freq == report_freq-1:
            edgeList.append(["report "+str(time-1)])
            print ("report:"+str(i)+'time:'+str(time))
            max_edge_count = max(max_edge_count,len(active_edge_set))

        if not dyn_input_param.incremental and i > sliding_window_size:
            delete_edge =edge_list[left_window][0]
            delete_edge.sort()
            if not dyn_input_param.duplicate:
                if str(delete_edge) in active_edge_set:
                    active_edge_set.pop(str(delete_edge))
                    # add the deletion symbol at the beginning
                    mod_edge = ['-']+delete_edge
                    edgeList.append(mod_edge)
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
            edgeList.append(mod_edge)

    edgeList.append(["report " + str(edge_list[-1][1])])
    n = len(vertex_set) # no of distinct vertices in the hypergraph
    print(max_edge_count,report_freq)
    if dbLPH.check_duplicate(edgeList):
        print("passed duplicate check")
    else:
        print("failed duplicate check")
    return edgeList,n,max_card,max_edge_count

def write_hg_file(edgeList,n,maxCard,max_edge_count,dyn_input_param):

    # print('n='+str(n)+'macCard='+str(maxCard))
    #write the edge list to a compressed text file
    with open (dyn_input_param.out_f_path, 'wt') as outFile:
        # Write header information
        outFile.write('% '+str(dyn_input_param)+'\n')
        outFile.write('% n = ' + str(n)+
                      ', max cardinality = '+str(maxCard)+
                      ',max_edge_count = '+str(max_edge_count)+'\n')
        for edge in edgeList:
            # Remove single author paper here
            str_edge = [str(element) for element in edge]
            outFile.write(' '.join(str_edge))
            outFile.write('\n')

def create_dynamic_hg(paper_list, dyn_input_param):
    edgeList,n,maxCard,max_edge_count = make_dynamic_hg(paper_list,dyn_input_param)

    # Now write the hg to a txt file
    write_hg_file(edgeList,n,maxCard,max_edge_count,dyn_input_param)

def read_edge_file(f_path):
    edge_list = []
    with open(f_path,'rt') as f:
        for line in f:
            if line[0]!= '%':
                edge = line.split()
                time = int(edge[-1])
                edge = [int(item) for item in edge[:-1]]
                edge_list.append([edge,time])
    return edge_list

if __name__== '__main__':

    dir_path = 'data/tags-stack-overflow/'
    input_file_path = dir_path + 'tags-stack-overflow-hg.txt'

    edgeList = read_edge_file(input_file_path)

    sliding_windows = [20]

    # Does the hypergraph support having multiple hyperedges?
    duplicate = False

    # If incremental is set to yes, then only hyperedge insertion will occur and
    # no hyperedge deletion will be done
    incremental = True

    for sliding_window_size in sliding_windows:
        #Output file names
        dyn = 'D'
        if incremental:
            dyn = 'I'

        hg_f_name = dir_path + '/tags_stack_of.hg.'+dyn+'.'+str(sliding_window_size)+'.txt'

        # Collect the input parameters for the dynamic graph
        dyn_input_param = DynInputParam(duplicate = duplicate,
                                            incremental = incremental,
                                            sliding_window_factor = sliding_window_size,
                                            out_f_path = hg_f_name
                                            )

        create_dynamic_hg(edgeList, dyn_input_param)

