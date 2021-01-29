from ortools.graph import pywrapgraph

def maxFlow (n,m,edgeList,edgeWeights,threshold):

def densestSubgraphFlow(n,m,edgeList,edgeWeights):
    lower = int(m/n)
    upper = sum(edgeWeights)
    threshold = 0
    optDensity = 0.0
    optS = []
    while upper - lower >= 1/n^2 :
        threshold = (upper+lower)*0.5
        flow_val,subsetS = maxFlow (n,m,edgeList,edgeWeights,threshold)
        if flow_val < m:
            lower = threshold
            #TODO store the optimal solution
        else:
            upper = threshold

def fileReader (hypergraph_file):

    batch = 1000000 # Read the first batch of lines from the file; for debug
    n = 0
    m = 0
    vertexSet = {}
    edgeSet = {}
    edgeList = []
    with open (hypergraph_file, 'rt') as file:
        line_no = 0
        for line in file:
            currentline = line.split(",")
            if line_no == 0:
                # n = int(currentline[0])
                # m = int(currentline[1])
                print(n,m)
            elif line_no < batch:
                if int(currentline[0]) <= 1: # ignore  self-loop for now TODO change?
                    continue
                edge = []
                for i in range(1,int(currentline[0])+1): # Read tge vertices in the hyperedge
                    u = int(currentline[i])
                    if u not in vertexSet:
                        vertexSet[u] = 1
                        n += 1
                    edge.append(u)
                if str(edge) not in edgeSet:
                    edgeList.append(edge)
                    edgeSet[str(edge)] = 1
                    # print(edge)
                    m += 1
            else:
                break
            line_no += 1
    edgeWeights = [1]* len(edgeList) # Considering unweighted case
    print(n,m,line_no,len(edgeList))
    return n,m,edgeList,edgeWeights

if __name__ == '__main__':
    hypergraph_filename = '../data/dblp.hypergraph.txt'
    n,m,edgeList,edgeWeights = fileReader(hypergraph_filename)
    densestSubgraphFlow(n,m,edgeList,edgeWeights)