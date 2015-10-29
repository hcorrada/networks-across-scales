from graph import Graph
from intlabels import IntLabel
from cycle_generator import SimpleGraphGenerator, CycleGenerator

g = Graph()
sources=range(10) + [2,6]
targets=[3,0,1,2,2,4,5,9,7,6,6,8]
for s,t in zip(sources,targets):
    source = g.get_or_make_node(IntLabel(s))
    target = g.get_or_make_node(IntLabel(t))
    g.add_edge(source, target)

print g.debug_print()

graph_generator = SimpleGraphGenerator(g)

node_to_bypass = graph_generator.find_bypass_node(g)
print "node to bypass"
print node_to_bypass.debug_print()
print

for new_graph in graph_generator.get_bypass_graphs(g, node_to_bypass):
    print "next graph"
    print new_graph.debug_print()
    print


graph_generator = SimpleGraphGenerator(g)
for graph in graph_generator:
    print "next graph:"
    print graph.debug_print()
    print

cycle_generator = CycleGenerator(g)
for cycle in cycle_generator:
    print "next cycle"
    print cycle
