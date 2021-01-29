from ortools.linear_solver import pywraplp
import os

# Parameters
# n: number of vertices
# m: number of edges
# edgeList: list of edges

def densestSubgraph(edgeList, edgeWeights):
    # Instantiate a Glop solver, naming it SolveDensest.
    solver = pywraplp.Solver('SolveDensest',
                             pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
    # solver.EnableOutput()
    # print(n,m)
    # Create the variables
    x = {}
    y = {}
    # Objective: maximize the weighted sum of edges \sum_{e} y_e w_e
    objective = solver.Objective()
    for i,edge in enumerate(edgeList):
        y[i] = solver.NumVar(0.0,solver.infinity(),'e'+str(i))
        objective.SetCoefficient(y[i], edgeWeights[i])
        for u in edge:
            if u not in x: # the variable corresponding to u is not yet created
                x[u] = solver.NumVar(0.0,solver.infinity(),'v'+str(u))
                objective.SetCoefficient(x[u], 0)
    print(len(x),len(y))
    objective.SetMaximization()

    # Create the constraints, r many per edge where r is the cardinality of the edge
    #  x_u - y_e >= 0 for each u \in e
    # One more constraint asserting that sum of x_u is 1
    constraints = []
    index = 0
    for i,edge in enumerate(edgeList):
        for u in edge:
            constraints.append(solver.Constraint(0, solver.infinity()))
            constraints[index].SetCoefficient(x[u],1)
            constraints[index].SetCoefficient(y[i],-1)
            index +=1
    constraints.append(solver.Constraint(0,1))
    for u in x:
        constraints[index].SetCoefficient(x[u],1)
    # print(index)
    status = solver.Solve()
    solution_set = []
    density_val = 0
    if status == solver.OPTIMAL:
        print('Optimal Found')
        print('Objective value =', solver.Objective().Value())
        # for i,edge in enumerate(edgeList):
        #     print(y[i].solution_value())
        subset_size = 0
        density_val = solver.Objective().Value()
        for u in x:
            if x[u].solution_value() > 0:
                subset_size += 1
                solution_set.append(u)
                # print(x[u].solution_value())
        print("Subset Size = ",subset_size)
    else:  # No optimal solution was found.
        if status == solver.FEASIBLE:
            print('A potentially suboptimal solution was found.')
        else:
            print('The solver could not solve the problem.')
    return solution_set, density_val

# The function to map the solution set of
# vertices to actual identifier for the vertices

def solutionMapper (solution_set, hypergraph_mapper_file):
    n = 0
    author_list = {}
    print(solution_set)
    with open(hypergraph_mapper_file, 'rt') as file:
        line_no = 0
        for line in file:
            line_no += 1
            if line_no == 1: # the first line constains total vertex number
                continue
            currentline = line.split(',')
            author_list[int(currentline[1])] = currentline[0]
    # print(author_list)
    # for x in solution_set:
    #     print(author_list[x])
    author_set = [author_list[x] for x in solution_set]
    author_set.sort()
    # print(author_set)


def fileReader (hg_file):
    n = 0
    m = 0
    vertexSet = {}
    edgeSet = {}
    edgeList = []
    with open (hg_file, 'rt') as file:
        line_no = 0
        for line in file:
            line_no += 1
            currentline = line.split(" ")
            if line_no == 1:
                print(currentline)
                continue
            if len(currentline) <= 2: # ignore  self-loop for now TODO change?
                continue
            edge = []
            # Read the vertices in the hyperedge
            # The last entry in each hyperedge is the year which we ignore
            for i in range(len(currentline)-1):
                u = int(currentline[i])
                if u not in vertexSet:
                    vertexSet[u] = 1
                    n += 1
                edge.append(u)
            # print(edge)
            if str(edge) not in edgeSet:
                edgeList.append(edge)
                edgeSet[str(edge)] = 1
                # print(edge)
                m += 1
    edgeWeights = [1]* len(edgeList) # Considering unweighted case
    print(n,m,line_no,len(edgeList))
    return edgeList,edgeWeights

def static_densest_subgraph(hg_static_file, hg_author_file):
    edgeList,edgeWeights = fileReader(hg_static_file)
    solution_set = densestSubgraph(edgeList,edgeWeights)
    solutionMapper (solution_set,hg_author_file)

# report the denest subgraph after every stat_interval
def dynmaic_densest_subgraph(hg_dyn_file, hg_author_file,out_f, out_dir):

    edgeSet = {}
    report_count = 0
    edgeList = []
    sol_set = []
    sol_val = []
    labels = []
    with open (hg_dyn_file, 'rt') as file:
        line_no = 0
        for line in file:
            line_no += 1
            currentline = line.split(" ")
            # ignore the lines that does not start with + or - or "report"
            if currentline[0] != '+' and currentline[0] != '-' and currentline[0]!= 'report':
                print(line_no,currentline)
                continue
            # if len(currentline) <= 2: # ignore  self-loop for now TODO change?
            #     print(line_no,currentline)
            #     continue
            # print(currentline[0])
            if currentline[0]== 'report':
                report_count += 1
                print("Running LP: "+str(report_count))
                edgeWeights = [1]* len(edgeList) # Considering unweighted cases
                solution_set, density_val = densestSubgraph(edgeList,edgeWeights)
                sol_set.append(solution_set)
                sol_val.append(density_val)
                labels.append(currentline[1].rstrip('\n'))
                # solutionMapper (solution_set,hg_author_file)
            else:
                # Read the vertices in the hyperedge
                # ignore the first character which is either a + or -
                edge = []
                for i in range(1,len(currentline)):
                    u = int(currentline[i])
                    edge.append(u)
                # The last entry in each hyperedge if it is
                # edge insertion is the year (which we ignore)
                if currentline[0] == '+':
                    edge = edge[:-1]
                    # if str(edge) not in edgeSet:
                    edgeList.append(edge)
                        # edgeSet[str(edge)] = 1
                elif currentline[0] == '-':
                    # if str(edge) in edgeSet:
                    if edge in edgeList: # Ensuring the hyperedge is present in the graph
                        edgeList.remove(edge)
                        # edgeSet.pop(str(edge))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(out_dir+out_f, 'wt') as out_file:
        out_file.write("Label Density Size \n")
        for i,label in enumerate(labels):
            out_file.write(label + " " + str(sol_val[i])
                           + "  " + str(len(sol_set[i]))+"\n")



if __name__ == '__main__':
    hg_filename = '../data/dblp.theory.static.hypergraph.txt'
    hg_author_filename = '../data/dblp.theory.static.hypergraph_author.txt'
    # static_densest_subgraph (hg_filename,hg_author_filename)

    # Input filename
    file_id = 'tags-math-sx'
    hg_dyn_f_list = [#'../data/threads_stack_overflow_full/th_stack_of.hg.D.-1.txt',
                     # '../data/tags-stack-overflow/tags_stack_of.hg.D.125.txt',
                     # '../data/DAWN/DAWN.hg.I.1.txt',
                     '../data/'+file_id+'/'+file_id+'.hg.I.1.txt',
                     # '../data/dblp/dblp.all.hg.D.1985.15.txt',
                     # '../data/dblp/dblp.all.hg.D.1985.15.txt'
                     ]
    hg_dyn_auth_f_list = ['../data/dblp.all.hg_author.dyn.1985.5.txt',
                          '../data/dblp.all.hg_author.dyn.1985.10.txt',
                          '../data/dblp.all.hg_author.dyn.1985.15.txt'
                          ]

    # Output directory and filename
    out_dir = '../data/output/'+file_id+'/'

    for i,hg_dyn_f in enumerate(hg_dyn_f_list):
        out_f = 'out.Exact.' + hg_dyn_f.rsplit('/')[-1]
        dynmaic_densest_subgraph(hg_dyn_f, hg_dyn_auth_f_list[i], out_f, out_dir)
