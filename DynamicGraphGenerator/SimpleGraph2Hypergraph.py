from typing import NamedTuple


# Input parameter class
class InputParam(NamedTuple):
    graph_f_path: str
    out_f_path: str
    self_loop: bool = False
    timestamp_scale: int = 1
    timestamp_threshold: int = 0

def raw_data_to_list(input_param: InputParam):
    edges = []
    vertex_map = {}
    index = 0
    max_cardinality = 0
    with open(input_param.graph_f_path) as f_path:
        for line in f_path:
            items = line.split()
            str_edge = items[:-1]
            edge = []
            for vertex in str_edge:
                vertex = int(vertex)
                if vertex not in vertex_map:
                    vertex_map[vertex] = index
                    index += 1
                edge.append(vertex_map[vertex])
            edge.sort()
            time = int(items[-1])
            time = int(time/input_param.timestamp_scale)  # Converts milliseconds into a quarter
            if input_param.self_loop or len(edge) > 1:
                if time > input_param.timestamp_threshold:
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

# A simulator to execute graph to hypergraph conversion in the desired format
# The file paths are specific to the graph files found at the
# https://snap.stanford.edu/data/c.html. (temporal networks)
# To use with other type of
# directory structure, modify accordingly.

if __name__ == '__main__':

    data_tag = 'sx-mathoverflow'

    dir_path = '../data/' + data_tag + '/'
    graph_filepath = dir_path + 'source/' + data_tag + '.txt'
    out_f_path = dir_path + data_tag + '-hg.txt'

    self_loop = False

    # millisecond to quarter: 7776000000=1000*3600*24*90
    # timescale = 7776000000
    # second to quarter: 3600*24*30*3
    timescale = 3600*24*30*3
    # timescale = 1 # Copy timestamp as it is from the input simplices

    # ignore all hyperedges before the timestamp_threshold
    timestamp_threshold = 1

    input_params = InputParam(graph_f_path=graph_filepath,
                              out_f_path=out_f_path,
                              self_loop=self_loop,
                              timestamp_scale=timescale,
                              timestamp_threshold=timestamp_threshold
                              )

    edgeList, n = raw_data_to_list(input_params)
