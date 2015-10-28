from graph import Graph
from intlabels import IntLabel
from cycle_generator import SimpleGraphGenerator

g = Graph()
sources=range(10) + [2,6]
targets=[3,0,1,2,2,4,5,9,7,6,6,8]
for s,t in zip(sources,targets):
    g.add_edge(IntLabel(s),IntLabel(t))

print g

graph_generator = SimpleGraphGenerator(g)

for graph in graph_generator:
    print graph
