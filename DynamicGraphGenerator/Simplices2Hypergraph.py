from pathlib import Path  # Only available in python>3.5
from itertools import islice
from typing import NamedTuple


# Input parameter class
class InputParam(NamedTuple):
    nvert_f_path: str
    simplices_f_path: str
    time_f_path: str
    out_f_path: str
    self_loop: bool = False
    timestamp_scale: int = 1


# Check if a directory exists, and create one if it does not
def create_directory(directory_path):
    Path(directory_path).mkdirs(parent=True, exist_ok=True)


# A function to convert a simplices file into a hypergraph
# The simplices files are taken from the following website:
# https://www.cs.cornell.edu/~arb/data/
# The output hypergraph file is in the following format
# %n=?,m=?, max_cardinality=?
# n1 n2 n3 timestamp
# .....
# .....
def raw_data_to_list(input_param: InputParam):
    edges = []
    vertex_map = {}
    index = 0
    max_cardinality = 0
    with open(input_param.nvert_f_path) as nvert_f, \
            open(input_param.simplices_f_path) as simplices_f, \
            open(input_param.time_f_path) as time_f:
        for line in nvert_f:
            edge = []
            count = int(line)
            vertices = islice(simplices_f, count)
            for vertex in vertices:
                vertex = int(vertex)
                if vertex not in vertex_map:
                    vertex_map[vertex] = index
                    index += 1
                edge.append(vertex_map[vertex])
            edge.sort()
            time = int(time_f.readline())
            time = int(time/input_param.timestamp_scale)  # Converts milliseconds into a quarter
            if input_param.self_loop or len(edge) > 1:
                edges.append([edge, time])
                max_cardinality = max(max_cardinality, len(edge))

    edges.sort(key=lambda k: k[1])

    with open(input_param.out_f_path, 'wt') as out_file:
        out_file.write('% n=' + str(index) + ' m=' + str(len(edges)) + ' max_cadinality='+str(max_cardinality)
                       + '\n')
        for edge in edges:
            str_edge = [str(element) for element in edge[0]]
            str_edge.append(str(edge[1]))
            out_file.write(' '.join(str_edge))
            out_file.write('\n')

    return edges, index


# A simulator to execute simplices to hypergraph conversion
# The file paths are specific to the simplices files found at the
# https://www.cs.cornell.edu/~arb/data/. To use with other type of
# directory structure, modify accordingly.
if __name__ == '__main__':

    data_tag = 'tags-math-sx'

    dir_path = '../data/' + data_tag + '/'
    nvert_filepath = dir_path + 'source/' + data_tag + '-nverts.txt'
    simplices_filepath = dir_path + 'source/' + data_tag + '-simplices.txt'
    time_filepath = dir_path + 'source/' + data_tag + '-times.txt'
    out_f_path = dir_path + data_tag + '-hg.txt'

    self_loop = False
    timescale = 7776000000  # millisecond to quarter: 7776000000=1000*3600*24*90

    input_params = InputParam(nvert_f_path=nvert_filepath,
                              simplices_f_path=simplices_filepath,
                              time_f_path=time_filepath,
                              out_f_path=out_f_path,
                              self_loop=self_loop,
                              timestamp_scale=timescale
                              )

    edgeList, n = raw_data_to_list(input_params)
