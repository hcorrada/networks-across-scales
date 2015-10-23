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
    current_item = cycle._head
    next_item = current_item._next_item

    while current_item._node != end and next_item._node != start:
        current_item = next_item
        next_item = current_item._next_item

    cycle._head = next_item
    cycle._tail = current_item
    return cycle

def find_eulerian_path(g, k, d):
    # find start and end nodes (also checks balance conditions)
    result = _find_start_end_nodes(g)
    if result is None:
        return None

    start, end = result

    # add edge to connect start and end node
    g.add_edge(end.label(), start.label())

    # get eulerian cycle
    cycle_finder = CycleFinder(g)
    cycle = cycle_finder.find_cycle()

    # remove extra edge
    path = get_path_from_cycle(cycle, start, end)

    # return path
    return path
