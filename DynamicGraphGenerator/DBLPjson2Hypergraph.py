import gzip, json
from typing import NamedTuple

class DynInputParam(NamedTuple):
    area: str
    conf: list
    duplicate: bool
    incremental: bool
    year_window: int
    start_year: int
    out_f_path: str
    out_f_auth_path: str


# Takes a json file with the following line format:
# ["conf/www/DemaineHMMRSZ14", ["Erik D. Demaine", "MohammadTaghi Hajiaghayi", "Hamid Mahini", "David L. Malec", "S. Raghavan", "Anshul Sawant", "Morteza Zadimoghaddam"], 2014],
# and outputs two .txt files.
# The first one is a hypergraph in the following format:
# n m # no of vertices and edge
# r,v1,v2,...,vr,year
# the lines are sorted by year
# The second one is a mapping between the vertex index and authors in the hypergraph

def json_to_list(json_gz_filename, years_threshold = 0, area = None, confs = None):
    lines = []
    # Read the json file and convert it into a sorted list
    with gzip.open (json_gz_filename, 'rt') as file:
        for line in file:
            if line.strip () in '[]': continue
            line = line.rstrip ().rstrip (',')

            # Extract the paper key, authors and the year from each line
            tag, authors, year = json.loads (line)

            # if either the year or authors entry is null, ignore this entry
            # Also, ignore single author paper
            if year is None or authors is None or len(authors)<=1 or year < years_threshold:
                # print(line)
                continue

            if area == 'all':# if no specific area is given, consider all papers
                lines.append([authors,year])
            else:# if specific area of research is given, then
                # parse the tag and look for matching conference name
                    for conf in confs:
                        if conf in tag and year >= years_threshold:
                            # print(tag, authors, year)
                            lines.append([authors,year])

    # Sort the paper entries in lines by year
    lines.sort(key=lambda k:k[1])
    return lines

def make_static_hg(paper_list):
    # paper_list is sorted by year
    # Hash the authors to integer values in the sorted list
    index = 0
    maxCard = 0
    edgeList = []
    author_map = {}

    for authors,year in paper_list:
        edge = []
        for author in authors:
            if author not in author_map:
                author_map[author] = index
                index +=1
            edge.append(author_map[author])

        maxCard = max(maxCard,len(edge))
        edge.append(year)
        edgeList.append(edge)
    n = index # no of distinct vertices in the hypergraph
    return edgeList,n,maxCard,author_map


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

    for authors,year in paper_list:
        if year not in edge_added_by_year:
            edge_added_by_year[year] = 0
        dup = False
        if year - start_year >= dyn_input_param.year_window:
            # Add a report statement
            edgeList.append(["report "+str(year-1)])
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
            mod_edge = ['+']+edge+[year]
            edgeList.append(mod_edge)

            edge_added_by_year[year] += 1
    edgeList.append(["report "+str(paper_list[-1][1])])
    n = index # no of distinct vertices in the hypergraph
    # print(edge_added_by_year)
    # print(edge_removed_by_year)
    if check_duplicate(edgeList):
        print("Passed sanity check")
    return edgeList,n,maxCard,author_map

    # A method to verify that there are no duplicate hyperedges
    # at any point in the edgeList
def check_duplicate(edgeList):
    active_edge_set = {}
    # with open('data/dblp/dblp.all.hg.D.1985.5.txt') as f:
    #     for line in f:
    #         edge = line.split()
    for edge in edgeList:
        if edge[0] == '+':
            n_edge = edge[1:-1]
            n_edge = [int(vertex) for vertex in n_edge]
            n_edge.sort()
            if str(n_edge) in active_edge_set:
                print ("Duplicate edge:" + str(edge))
                return False
            else:
                active_edge_set[str(n_edge)] = 1
        elif edge[0] == '-':
            n_edge = edge[1:]
            n_edge = [int(vertex) for vertex in n_edge]
            n_edge.sort()
            if str(n_edge) not in active_edge_set:
                print ("Removing non existent edge:" + str(edge))
                return False
            else:
                del active_edge_set[str(n_edge)]
    return True

def write_hg_file(edgeList,n,maxCard,dyn_input_param):

    # print('n='+str(n)+'macCard='+str(maxCard))
    #write the edge list to a compressed text file
    with open (dyn_input_param.out_f_path, 'wt') as outFile:
        # Write header information
        outFile.write('% '+str(dyn_input_param)+'\n')
        outFile.write('% n = ' + str(n)+', max cardinality = '+str(maxCard)+'\n')
        for edge in edgeList:
            # Remove single author paper here
            str_edge = [str(element) for element in edge]
            outFile.write(' '.join(str_edge))
            outFile.write('\n')

def write_author_map (hg_author_filename,authorList,n):
    with open(hg_author_filename, 'wt') as authorFile:
        authorFile.write('Number of Authors = '+ str(n)+'\n')
        for author in authorList:
            authorFile.write(author +',' + str(authorList[author]) + '\n')

def create_dynamic_hg(paper_list, dyn_input_param):
    edgeList,n,maxCard,author_map = make_dynamic_hg(paper_list,dyn_input_param)

    # Now write the hg to a txt file
    write_hg_file(edgeList,n,maxCard,dyn_input_param)

    # write the author name to index mapping in a txt file
    write_author_map(dyn_input_param.out_f_auth_path,author_map,n)


def create_static_hg(paper_list,hg_f_name, hg_auth_f_name):
    edgeList,n,maxCard,author_map = make_static_hg(paper_list)

    # Now write the hg to a txt file
    write_hg_file(edgeList,n,maxCard,hg_f_name)

    # write the author name to index mapping in a txt file
    write_author_map(hg_auth_f_name,author_map,n)



if __name__ == '__main__':
    json_gz_filename = 'data/dblp.conf_journal.json.gz'

    years_threshold = 1985

    # area = 'theory'
    # confs = ['/soda/','/stoc/','/focs/','/icalp/', '/coco/', #coco== CCC
    #          '/podc/','/esa/','/mfcs',
    #          '/stacs/','/ipco/','/pods/']

    area = 'all'
    confs = []

    # year_window dictates the no of years in consideration for the hypergraph
    # The desnity will be reported for each sliding window of year_window many years
    # For eg, consider year_window = 5 and incremental = false
    # first the papers from 1985-1989 will be added. There will a report command
    # Then, papers from 1990 will be added
    # and 1985 will be removed.
    year_windows = [5,10,15]

    # Does the hypergraph support having multiple hyperedges?
    duplicate = False

    # If incremental is set to yes, then only hyperedge insertion will occur and
    # no hyperedge deletion will be done
    incremental = False

    # Parse the json file and create a list where each
    # element corresponds to a paper: it contains the author list and the year
    # The list is sorted by year
    # Further, the list only contains paper from the confs in areas which
    # appear after the years-threshold
    paper_list = json_to_list(json_gz_filename, years_threshold, area, confs)

    for year_window in year_windows:
        #Output file names
        dyn = 'D'
        if incremental:
            dyn = 'I'

        hg_f_name = 'data/dblp/dblp.'+area+'.hg.'+dyn+'.'+str(years_threshold)+'.'+str(year_window)+'.txt'
        hg_auth_f_name = 'data/dblp/dblp.'+area+'.hg_author.'+dyn+'.'+str(years_threshold)+'.'+str(year_window)+'.txt'

        # Collect the input parameters for the dynamic graph
        dyn_input_param = DynInputParam(area = area, conf = confs,
                                        duplicate = duplicate,
                                        incremental = incremental,
                                        year_window= year_window,
                                        start_year = years_threshold,
                                        out_f_path = hg_f_name,
                                        out_f_auth_path = hg_auth_f_name)

        create_dynamic_hg(paper_list, dyn_input_param)


    # if check_duplicate([]):
    #     print("Passed Sanity Check")


    # Create a static hypegraph
    # create_static_hg(paper_list,hg_f_name, hg_auth_f_name)
