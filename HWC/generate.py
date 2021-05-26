import random

filename = "test.in"
def random_graph(seed = 0, n = 3000,
                dl = 2,
                dr = 2,
                opNumber = 10000,
                addProb = 0.6,
                i = 0,
                valid = set(),
                edges = [],
                edge_id = 0):
    random.seed(a=seed);
    total_edges = set()
    with open(filename, 'w') as f:
        vertices = range(0, n)
        while i < opNumber:
            #print i
            l = []
            if random.random() < addProb:
                l.append('+')
                tried = 1
                while True: 
                    edge = random.sample(vertices, random.randint(dl, dr))
                    edge.sort()
                    te = tuple(edge)
                    if te not in total_edges:
                        total_edges.add(te)
                        break
                    tried += 1
                    if tried > 30:
                        vertices = range(1, n)
                edges.append(edge)
                l.extend(edge)
                valid.add(edge_id)
                edge_id += 1
            else:
                l.append('-')
                if len(valid) == 0:
                    continue
                j = random.sample(valid, 1)[0]
                valid.remove(j)
                l += edges[j]
                # l += [j]
            print >> f, " ".join(map(str, l))
            i += 1
            if i % 40 == 0:
                start = random.randint(0, n - 100)
                vertices = range(start, start+random.randint(20, 30))

# seed [47]
if __name__ == '__main__':
    random_graph(seed = 47)
