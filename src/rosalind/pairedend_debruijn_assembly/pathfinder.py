class DummyPath:
    def get_string(self):
        return None

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

def find_eulerian_path(g):
    # find start and end nodes (also checks balance conditions)
    result = _find_start_end_nodes(g)
    if result is None:
        return None

    start, end = result

    # add edge to connect start and end node

    # get eulerian cycle

    # remove extra edge

    # return path
    return DummyPath()
