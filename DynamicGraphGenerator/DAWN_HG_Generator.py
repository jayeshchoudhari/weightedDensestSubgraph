from typing import NamedTuple
from DynamicGraphGenerator import DBLPjson2Hypergraph as dbLPH


class DynInputParam(NamedTuple):
    duplicate: bool
    incremental: bool
    year_window: int
    out_f_path: str

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

def make_dynamic_hg(paper_list, dyn_input_param):


    # paper_list is sorted by year
    # Hash the authors to integer values in the sorted list
    index = 0
    maxCard = 0
    edgeList = []
    author_map = {}
    active_edge_set = {}
    paper_by_year = {}
    start_year = paper_list[0][1]

    edge_added_by_year = {}
    edge_removed_by_year = {}

    max_card = 0
    max_edge_count = 0

    for authors,year in paper_list:
        if year not in edge_added_by_year:
            edge_added_by_year[year] = 0
        dup = False
        if year - start_year >= dyn_input_param.year_window:
            # Add a report statement
            edgeList.append(["report "+str(year-1)])
            max_edge_count = max(max_edge_count,len(active_edge_set))
            if year not in edge_removed_by_year:
                edge_removed_by_year[year] = 0
            # For fully-dynamic input graph, remove edges
            if not dyn_input_param.incremental:
                for edge in paper_by_year[start_year]:
                    if not dyn_input_param.duplicate:
                        active_edge_set.pop(str(edge))
                    # add the deletion symbol at the beginning
                    # edge.insert(0,'-')
                    mod_edge = ['-']+edge
                    edgeList.append(mod_edge)
                    edge_removed_by_year[year] += 1

            start_year += 1

        edge = []
        for author in authors:
            if author not in author_map:
                author_map[author] = index
                index +=1
            edge.append(author_map[author])
            edge.sort()

        if str(edge) not in active_edge_set:
            active_edge_set[str(edge)] = 1
        else:
            dup = True

        if dyn_input_param.duplicate or not dup:
            maxCard = max(maxCard,len(edge))
            if year not in paper_by_year:
                paper_by_year[year] = [edge]
            else:
                paper_by_year[year].append(edge)
            # edge.insert(0,'+')
            # edge.append(year)
            max_card = max(max_card,len(edge))
            mod_edge = ['+']+edge+[year]
            edgeList.append(mod_edge)

            edge_added_by_year[year] += 1
    edgeList.append(["report "+str(paper_list[-1][1])])
    n = index # no of distinct vertices in the hypergraph
    # print(edge_added_by_year)
    # print(edge_removed_by_year)
    if dbLPH.check_duplicate(edgeList):
        print("Passed sanity check")
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


def create_dynamic_hg(edge_list, dyn_input_param):
    edgeList,n,maxCard,max_edge_count = make_dynamic_hg(edge_list,dyn_input_param)

    # Now write the hg to a txt file
    write_hg_file(edgeList,n,maxCard,max_edge_count,dyn_input_param)


if __name__ == '__main__':


    file_id = 'tags-math-sx'

    dir_path = 'data/'+file_id+'/'
    input_file_path = dir_path + file_id + '-hg.txt'

    edgeList = read_edge_file(input_file_path)


    # year_window dictates the no of years in consideration for the hypergraph
    # The desnity will be reported for each sliding window of year_window many years
    # For eg, consider year_window = 5 and incremental = false
    # first the papers from 1985-1989 will be added. There will a report command
    # Then, papers from 1990 will be added
    # and 1985 will be removed.
    year_windows = [1]

    # Does the hypergraph support having multiple hyperedges?
    duplicate = False

    # If incremental is set to yes, then only hyperedge insertion will occur and
    # no hyperedge deletion will be done
    incremental = False

    for year_window in year_windows:
        #Output file names
        dyn = 'D'
        if incremental:
            dyn = 'I'

        hg_f_name = dir_path + file_id + '.hg.'+dyn+'.'+str(year_window)+'.txt'


        # Collect the input parameters for the dynamic graph
        dyn_input_param = DynInputParam(duplicate = duplicate,
                                        incremental = incremental,
                                        year_window= year_window,
                                        out_f_path = hg_f_name)

        create_dynamic_hg(edgeList, dyn_input_param)
