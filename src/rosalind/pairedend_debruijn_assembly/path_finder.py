from cycle_finder import CycleFinder

def _find_start_end_nodes(g):
    # compute node degrees
    node_degrees = g.node_degrees()

    # check balance conditions
    start = end = None
    for node in g:
        in_degree, out_degree = node_degrees[node.label()]

        if in_degree == out_degree:
            continue

        if in_degree == out_degree + 1 and end is None:
            end = node
            continue

        if in_degree == out_degree - 1 and start is None:
            start = node
            continue

        return None
    if start is None or end is None:
        return None
    return start, end

def get_path_from_cycle(cycle, start, end):
    end_items = cycle.find_all_items(end)

    start_item = None
    end_item = None

    for item in end_items:
        end_item = item
        start_item = item._next_item

        if start_item._node == start:
            break
        else:
            start_item = None
    assert(start_item is not None)
    assert(end_item is not None)

    cycle._tail._next_item = cycle._head._next_item

    cycle._head = start_item
    cycle._tail = end_item
    cycle._head._prev_item = None
    cycle._tail._next_item = None

    return cycle

def find_eulerian_path(g):
    # find start and end nodes (also checks balance conditions)
    result = _find_start_end_nodes(g)
    if result is None:
        return None

    start, end = result

    # add edge to connect start and end node
    g.add_edge(end.label(), start.label())

    # get eulerian cycle
    cycle_finder = CycleFinder(g)
    cycle = cycle_finder.run()

    # remove extra edge
    path = get_path_from_cycle(cycle, start, end)

    # return path
    return path
