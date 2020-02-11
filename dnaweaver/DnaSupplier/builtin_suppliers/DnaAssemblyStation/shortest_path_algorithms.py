"""Optimization techniques"""

import numpy as np
import networkx as nx
from heapq import heappush, heappop
from itertools import count


class NoSolutionFoundError(Exception):
    """Error thrown when an optimization problem seems to have no solution."""
    pass


def astar_path(G, source, target, heuristic=None, weight=None):
    """Return a list of nodes, A* shortest path between source and target.

    Uses the A* ("A-star") algorithm.

    There may be more than one shortest path.  This returns only one.

    This function is taken from the Networkx project, with modifications for
    lazy weight computations (weight can now be a function).
    This function is therefore placed under the BSD licence.

    TO DO: submit that change to the NetworkX project so that this rewrite
    can be avoided.


    Parameters
    ----------

    G : NetworkX graph
    source : node
       Starting node for path
    target : node
       Ending node for path
    heuristic : function
       A function to evaluate the estimate of the distance
       from the a node to the target.  The function takes
       two nodes arguments and must return a number.
    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight.

    Raises
    ------

    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> print(nx.astar_path(G,0,4))
    [0, 1, 2, 3, 4]
    >>> G=nx.grid_graph(dim=[3,3])  # nodes are two-tuples (x,y)
    >>> def dist(a, b):
    ...    (x1, y1) = a
    ...    (x2, y2) = b
    ...    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    >>> print(nx.astar_path(G,(0,0),(2,2),dist))
    [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)]

    See Also
    --------
    shortest_path, dijkstra_path
    """
    if source not in G or target not in G:
        msg = 'Either source {} or target {} is not in G'
        raise nx.exception.NetworkXNoPath(msg.format(source, target))

    if heuristic is None:
        # The default heuristic is h=0 - same as Dijkstra's algorithm
        def heuristic(u, v):
            return 0
    if weight is None:
        def weight(*a):
            return 1
    elif isinstance(weight, str):
        key = weight
        def weight(n1, n2, props):
            return props[key]

    push = heappush
    pop = heappop
    c = count()
    queue = [(0, next(c), source, 0, None)]
    enqueued = {}
    explored = {}

    while queue:
        _, __, curnode, dist, parent = pop(queue)

        if curnode == target:
            path = [curnode]
            node = parent
            while node is not None:
                path.append(node)
                node = explored[node]
            path.reverse()
            return path

        if curnode in explored:
            continue

        explored[curnode] = parent

        for neighbor, w in G[curnode].items():
            if neighbor in explored:
                continue
            ncost = dist + weight(curnode, neighbor, w)
            if ncost == np.inf:
                continue
            if neighbor in enqueued:
                qcost, h = enqueued[neighbor]
                if qcost <= ncost:
                    continue
            else:
                h = heuristic(neighbor, target)
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))

    raise nx.NetworkXNoPath("Node %s not reachable from %s" % (source, target))


def shortest_valid_path(graph, start, end, nodes_constraints=(),
                        size_limit=None, min_step=None,
                        compatibility_search_cutoff=100,
                        penalty=0):
    if size_limit is not None:
        def f(penalty):
            try:
                path = shortest_valid_path(
                    graph, start, end, penalty=penalty, size_limit=None,
                    nodes_constraints=nodes_constraints,
                    compatibility_search_cutoff=100
                )
                return path, len(path) - 1
            except NoSolutionFoundError:
                return None, np.inf

        weights = [data["weight"] for _, _, data in graph.edges(data=True)]
        if weights == []:
            raise NoSolutionFoundError()
        penalty = max(weights)
        step = -0.5 * penalty

        best_path, best_path_size = f(penalty)
        if best_path_size > size_limit:
            raise NoSolutionFoundError()

        while abs(step) > min_step:
            penalty = penalty + step
            path, path_size = f(penalty)
            if path_size == size_limit:
                return path
            elif path_size > size_limit:
                step = 0.5 * abs(step)
            else:  # path_size < size_limit
                step = -0.5 * abs(step)
                best_path = path
        return best_path

    elif penalty > 0:
        _, _, edge_datas = zip(*graph.edges(data=True))
        for data in edge_datas:
            data["weight"] += penalty
        try:
            path = shortest_valid_path(
                graph, start, end, penalty=0, size_limit=None,
                compatibility_search_cutoff=compatibility_search_cutoff,
                nodes_constraints=nodes_constraints
            )
            for data in edge_datas:
                data["weight"] -= penalty
        except nx.NetworkXNoPath:
            for data in edge_datas:
                data["weight"] -= penalty
            raise NoSolutionFoundError()
        return path

    elif len(nodes_constraints) == 0:
        return nx.dijkstra_path(graph, start, end, weight='weight')

    else:
        def path_is_valid(path):
            return all(nodes_constraint(path)
                       for nodes_constraint in nodes_constraints)
        tentative_shortest_path = nx.dijkstra_path(graph, start, end,
                                                   weight='weight')
        if path_is_valid(tentative_shortest_path):
            return tentative_shortest_path
        shortest_paths = nx.shortest_simple_paths(graph, start, end,
                                                  weight='weight')
        for i in range(compatibility_search_cutoff):
            shortest_path = next(shortest_paths)
            if path_is_valid(shortest_path):
                return shortest_path
        raise nx.NetworkXNoPath("Could not find a solution verifying the cuts"
                                 " set constraints, after %d tries." %
                                 compatibility_search_cutoff)
