from graph import Graph
from intlabels import IntLabel
from path_finder import find_eulerian_path

g = Graph()
sources=[0,1,2,3,3,6,6,7,8,9]
targets=[2,3,1,0,4,3,7,8,9,6]

for s,t in zip(sources,targets):
    source = g.get_or_make_node(IntLabel(s))
    target = g.get_or_make_node(IntLabel(t))
    g.add_edge(source, target)

print g

path = find_eulerian_path(g)

print path
