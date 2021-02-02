import gzip
import json
from typing import NamedTuple


# Input parameter class
class InputParam(NamedTuple):
    json_f_path: str
    area: str
    confs: list
    years_threshold: int
    out_f_path: str
    out_author_f_path: str


# Takes a json file with the following line format:
# ["conf/www/DemaineHMMRSZ14", ["Erik D. Demaine", "MohammadTaghi Hajiaghayi", "Hamid Mahini", "David L. Malec", "S. Raghavan", "Anshul Sawant", "Morteza Zadimoghaddam"], 2014],
# and outputs two .txt files.
# The first one is a hypergraph in the following format:
# n m # no of vertices and edge
# r,v1,v2,...,vr,year
# the lines are sorted by year
# The second one is a mapping between the vertex index and authors in the hypergraph
def json_to_list(input_param: InputParam):

    lines = []
    vertex_map = {}
    index = 0
    max_cardinality = 0

    # Read the json file and convert it into a sorted list
    with gzip.open(input_param.json_f_path, 'rt') as file:
        for line in file:
            if line.strip() in '[]': continue
            line = line.rstrip().rstrip(',')

            # Extract the paper key, authors and the year from each line
            tag, authors, year = json.loads (line)

            # if either the year or authors entry is null, ignore this entry
            # Also, ignore single author paper
            if year is None or authors is None or len(authors)<=1 or year < input_param.years_threshold:
                # print(line)
                continue
            edge = []
            for author in authors:
                if author not in vertex_map:
                    vertex_map[author] = index
                    index += 1
                edge.append(vertex_map[author])
            edge.sort()

            if input_param.area == 'all':# if no specific area is given, consider all papers
                lines.append([edge, year])
                max_cardinality = max(max_cardinality, len(edge))
            else:   # if specific area of research is given, then
                    # parse the tag and look for matching conference name
                for conf in input_param.confs:
                    if conf in tag and year >= input_param.years_threshold:
                        # print(tag, authors, year)
                        lines.append([edge, year])
                        max_cardinality = max(max_cardinality, len(edge))

    # Sort the paper entries in lines by year
    lines.sort(key=lambda k: k[1])

    with open(input_param.out_f_path, 'wt') as out_file:
        out_file.write('% n=' + str(index) + ' m=' + str(len(lines)) + ' max_cadinality='+str(max_cardinality)
                       + '\n')
        for edge in lines:
            str_edge = [str(element) for element in edge[0]]
            str_edge.append(str(edge[1]))
            out_file.write(' '.join(str_edge))
            out_file.write('\n')

    with open(input_param.out_author_f_path, 'wt') as authorFile:
        authorFile.write('Number of Authors = ' + str(index)+'\n')
        for author in vertex_map:
            authorFile.write(str(vertex_map[author]) + ',' + author + '\n')

    return lines, index


if __name__ == '__main__':

    json_gz_filename = '../data/dblp.conf_journal.json.gz'
    dir_path = '../data/dblp/'

    # area = 'theory'
    # confs = ['/soda/','/stoc/','/focs/','/icalp/', '/coco/', #coco== CCC
    #          '/podc/','/esa/','/mfcs',
    #          '/stacs/','/ipco/','/pods/']

    area = 'all'
    confs = []

    years_threshold = 1985
    out_f_path = dir_path + 'dblp.' + area + '.' + str(years_threshold) + '-hg.txt'
    out_auth_f_name = dir_path + 'dblp.author.' + area + '.' + str(years_threshold) + '-hg.txt'

    # Collect the input parameters for the hypergraph
    input_params = InputParam(json_f_path=json_gz_filename,
                              area=area,
                              confs=confs,
                              years_threshold=years_threshold,
                              out_f_path=out_f_path,
                              out_author_f_path=out_auth_f_name,
                              )
    # Parse the json file and create a list where each
    # element corresponds to a paper: it contains the author list and the year
    # The list is sorted by year
    # Further, the list only contains paper from the confs in areas which
    # appear after the years-threshold
    edge_list, n = json_to_list(input_params)
