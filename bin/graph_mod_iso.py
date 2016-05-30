#!/usr/bin/env python
import networkx as nx
import community as com
import numpy as np
import sys
import os

if len(sys.argv) != 2:
    print 'Usage: [input graphml]'
    sys.exit(1)

g = nx.read_graphml(sys.argv[1])

# no self-loops
g.remove_edges_from(g.selfloop_edges())

if g.order() == 0:
    print 'The graph contains no nodes'
    sys.exit(1)
else:
    try:
        # determine Louvain modularity score of entire graph
        if g.size() == 0:
            # no edges, modularity is 1
            mod = 1.0
        else:
            part = com.best_partition(g)
            mod = com.modularity(part, g)
    except Exception as e:
        print 'An exception occurred during modularity analysis'
        exc_type, exc_obj, exc_tb = sys.exc_info()        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]        print exc_type, fname, exc_tb.tb_lineno
        sys.exit(1)

# count isolates
n_iso = len(nx.isolates(g))

# mean degree
mean_deg = np.mean(nx.degree(g).values())

# median degree
med_deg = np.median(nx.degree(g).values())

print '#file #nodes #iso #inclu Modularity MeanDeg MedianDeg'
print '{0} {1} {2} {3} {4} {5} {6}'.format(sys.argv[1], g.order(), n_iso, g.order() - n_iso, mod, mean_deg, med_deg)
