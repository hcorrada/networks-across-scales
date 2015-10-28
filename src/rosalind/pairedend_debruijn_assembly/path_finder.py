from cycle_finder import CycleFinder
from cycle_generator import CycleGenerator

# find nodes corresponding to start and end of the path
# also checks for balancing conditions to find
# an Eulerian path
#
# the start node must have in-degree = out-degree - 1
# the end node must have in-degree = out-degree + 1
# all other nodes must have in-degree = out-degree
#
# input:
#   g: object of class Graph
# returns:
#   tuple of nodes: start, end if conditions are satisfied
#   otherwise returns None
def _find_start_end_nodes(g):
    # compute node degrees
    node_degrees = g.node_degrees()

    # check balance conditions
    start = end = None
    for node in g:
        in_degree, out_degree = node_degrees[node.label()]

        # node is balanced
        if in_degree == out_degree:
            continue

        # this is the first end node we find
        if in_degree == out_degree + 1 and end is None:
            end = node
            continue

        # this is the first start node we find
        if in_degree == out_degree - 1 and start is None:
            start = node
            continue

        # node is unbalanced, doesn't qualify for start or end
        # or we already have start or end node
        # conditions not met
        return None
    # check we found both start and end nodes
    if start is None or end is None:
        return None
    return start, end

# remove added edge from end to start node
# to recover Eulerian path from Eulerian cycle
#
# input:
#   cycle: object of class Path
#   start: Node object for start node in path
#   end: Node object for end node in path
#
# returns:
#   object of class Path
def get_path_from_cycle(cycle, start, end):
    # the idea is to find the edge in the path
    # connecting end -> start

    # find all items in the path that correspond
    # to the end node
    end_items = cycle.find_all_items(end)

    # now let's find which of those connects
    # to the start node
    start_item = None
    end_item = None

    for item in end_items:
        end_item = item
        start_item = item._next_item

        if start_item._node == start:
            break
        else:
            start_item = None

    # make sure we did the right thing
    assert(start_item is not None)
    assert(end_item is not None)

    # first of all, this is a cycle
    # so the head and tail point to the same node
    # so let's splice out the copy pointed to by the
    # head pointer in the cycle
    cycle._tail._next_item = cycle._head._next_item

    # now start the path to return at the start item
    cycle._head = start_item

    # end the path at the end item
    cycle._tail = end_item

    # and make sure head and tail are valid
    cycle._head._prev_item = None
    cycle._tail._next_item = None

    return cycle

# preprocess graph before findign Eulerian cycle
#
# 1. checks balance conditions
# 2. finds start and end nodes
# 3. adds end->start edge
#
# input:
#   g: object of class Graph
# returns:
#   tuple if graph satisfies conditions: new_graph, start_node, end_node
#   None otherwise
def _preprocess_graph(g):
    # find start and end nodes
    # if balance conditions for graph are not satisfied
    # this will return None
    result = _find_start_end_nodes(g)

    # check if balancing conditions are met
    if result is None:
        return None

    # unpack start and end nodes
    start, end = result

    # add edge to connect start and end nodes
    g.add_edge(end.label(), start.label())
    return g, start, end

# finds an Eulerian path, if it exists, in the given
# graph
#
# input:
#   g: object of class Graph
# returns:
#   object of class Path
def find_eulerian_path(g):
    # preprocess graph
    g, start, end = _preprocess_graph(g)

    # find Eulerian cycle
    cycle_finder = CycleFinder(g)
    cycle = cycle_finder.run()

    # remove extra edge from cycle to
    # recover Eulerican path
    path = get_path_from_cycle(cycle, start, end)

    # return path
    return path

# find valid eulerian path given graph
#
# input:
#   g: object of class Graph
#   validity_func: function Path->boolean that checks if generated path is valid
def find_valid_eulerian_path(g, validity_func):
    g, start, end = _preprocess_graph(g)

    # generate all Eulerian cycles
    cycle_generator = CycleGenerator(g)
    for cycle in cycle_generator:
        # get the next candidate cycle, turn into path
        path = get_path_from_cycle(cycle, start, end)
        if validity_func(path):
            # this is a valid path, so return it
            return path
    # couldn't find a valid path, so return None
    return None
