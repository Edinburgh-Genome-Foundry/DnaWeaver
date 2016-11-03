"""Optimization techniques"""

import numpy as np
from tqdm import tqdm
import itertools as itt
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
        weight = lambda *a: 0

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
            if neighbor in enqueued:
                qcost, h = enqueued[neighbor]
                if qcost <= ncost:
                    continue
            else:
                h = heuristic(neighbor, target)
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))

    raise nx.NetworkXNoPath("Node %s not reachable from %s" % (source, target))


def shortest_compatible_path(graph, start, end, nodes_constraints=(),
                             compatibility_search_cutoff=100):

    if len(nodes_constraints) == 0:
        return nx.dijkstra_path(graph, start, end)
    shortest_paths = nx.shortest_simple_paths(graph, start, end)
    for i in range(compatibility_search_cutoff):
        shortest_path = next(shortest_paths)
        if all([nodes_constraint(shortest_path)
                for nodes_constraint in nodes_constraints]):
            return shortest_path
    return nx.NetworkXNoPath("Could not find a solution verifying the cuts"
                             " set constraints, after %d tries." %
                             compatibility_search_cutoff)


def shortest_compatible_path_with_penalty(graph, start, end, penalty=0,
                                          nodes_constraints=(),
                                          compatibility_search_cutoff=100):
    """Compute the graph's shortest path after adding a penalty to each edge.

    A copy of the graph is done and each edge weight receives a additional
    constant penalty. As a result, paths with less edges are more likely to
    be the shortest paths with respect to total weight.

    Parameters:
    -----------

    graph
      A networkX directed graph

    start
      Starting node for the path

    end
      Target node for the path

    penalty
      Additional weight applied to all edges before computing the shortest
      path. The original graph remains unchanged.

    Returns:
    --------

    A path [node1, node2, node3...], or raises NoSolutionFoundError if no
    path was found

    """
    _, _, edge_datas = zip(*graph.edges(data=True))
    for data in edge_datas:
        data["weight"] += penalty
    try:
        path = shortest_compatible_path(
            graph, start, end,
            compatibility_search_cutoff=compatibility_search_cutoff,
            nodes_constraints=nodes_constraints
        )
    except nx.NetworkXNoPath:
        for data in edge_datas:
            data["weight"] -= penalty
        raise NoSolutionFoundError()
    return path


def shortest_compatible_path_with_size_limit(graph, start, end, size_limit,
                                             min_step, nodes_constraints=(),
                                             compatibility_search_cutoff=100):
    """Compute the graph's shortest path whose number of edges is below a limit

    The algorithm applies various penalties to the graph's weights to reduce
    the number of edges in the shortest path. It uses a bisection search to
    find the optimal penalty (yielding a shortest path with a number of edges
    smaller but as close as possible from `size_limit`)

    Parameters:
    -----------

    graph
      A networkX directed graph

    start
      Starting node for the path

    end
      Target node for the path

    size_limit
      Maximal number of edges in acceptable paths

    min_step
      Step size after which to stop the penalty bisection.

    Returns:
    --------

    A path [node1, node2, node3...], or raises NoSolutionFoundError if no
    path was found


    """
    def f(penalty):
        path = shortest_compatible_path_with_penalty(
            graph, start, end, penalty,
            nodes_constraints=nodes_constraints,
            compatibility_search_cutoff=100
        )
        if path == False:  # means no solution found
            return False
        return path, len(path) - 1

    penalty = max(data["weight"] for _, _, data in graph.edges(data=True))
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


def optimize_cuts_with_graph(sequence_length, segment_score_function,
                             cut_location_constraints=(),
                             segment_constraints=(), forced_cuts=(),
                             suggested_cuts=(),
                             cuts_set_constraints=(),
                             max_segment_length=None, a_star_factor=0,
                             path_size_limit=None, path_size_min_step=0.001,
                             progress_bars=False):
    """Find the sequence cuts which optimize the sum of segments scores.

    This is a very generic method meant to be applied to any sequence cutting
    optimization problem.

    Parameters
    ----------

    sequence_length
      Length of the sequence

    segment_score_function
      A function f( (start, end) ) -> score where (start, end) refers to the
      sequence segment start:end, and score is a float. The algorithm will
      produce cuts which minimize the total scores of the segments.

    cut_location_constraints
      List or tuple of functions `int -> bool` which return for each location
      whether the location should be considered as a cutting site (True) or
      discarded (False) in the decomposition of the sequence.
      The locations considered are the locations which pass every filter
      in the `cut_location_constraints` list.

    segment_constraints
      List or tuple of functions `(int, int) -> bool` which return for each
      segment `(start, end)` whether the segment is a valid segment for the
      decomposition of the sequence or whether it should be forbidden.
      The segmentss considered are the segments which pass every filter
      in the `segments_filters` list.

    forced_cuts
      List of locations at which the decomposition must imperatively cut, even
      if these cuts do not comply with the `cut_location_constraints`.

    max_segment_length
      Maximal length of the segments. Even though this could be specified using
      a filter in `segments_filters`, in practice providing this parameter
      accelerates computations as it allows to reduce the number
      of segments considered.

    a_star_factor
      If 0, the classical Dijkstra algorithm is used for path finding. Else,
      the a_star algorithm is used with a heuristic `h(x)=a_start_factor*(L-x)`
      where x is the location of a cutting point and L the length of the
      sequence. See the original DNAWeaver article for clearer explanations.
      Using a high A* factor can improve computing times several folds but
      yields suboptimal decompositions.

    path_size_limit
      Maximal number of edges for acceptable paths. None means no limit.
      Only works with `a_star_factor=0`

    path_size_min_step
      Minimal step for the bisection search when `path_size_limit` is not None
      (see `dijkstra_path_with_size_limit`)


    Returns
    -------

    graph
      The graph of the problem (a Networkx DiGraph whose nodes are indices of
      locations in the sequence and the "weight" of edge [i][j] is given by
      segment_score_function(segment[i][j])).

    best_cuts
      The list of optimal cuts (includes 0 and len(sequence)), which minimizes
      the total score of all the segments corresponding to the cuts.

    """

    forced_cuts = list(forced_cuts)
    suggested_cuts = list(suggested_cuts)

    if max_segment_length is None:
        max_segment_length = sequence_length

    # List all nodes, add 0 and seqlength, add forced+suggested cuts
    nodes = sorted(list(set([0, sequence_length] + [
        node
        for node in range(0, sequence_length)
        if all(fl(node) for fl in cut_location_constraints)
    ] + forced_cuts + suggested_cuts)))

    if forced_cuts != []:
        def forced_cuts_filter(segment):
            start, end = segment
            return not any((start < cut < end) for cut in forced_cuts)
        segment_constraints = [forced_cuts_filter] + list(segment_constraints)

    segments = []
    for i, start in enumerate(nodes if not progress_bars else
                              tqdm(nodes, desc="Filtering edges")):
        for end in nodes[i + 1:]:
            if end - start > max_segment_length:
                break
            elif all(fl(start, end) for fl in segment_constraints):
                segments.append((start, end))
    graph = nx.DiGraph()

    if (len(cuts_set_constraints) > 0) or (a_star_factor == 0):

        for start, end in (segments if not progress_bars else
                           tqdm(segments, desc='Computing edges')):
            weight = segment_score_function((start, end))
            if weight >= 0:
                graph.add_edge(start, end, weight=weight)

        try:
            if path_size_limit is not None:

                best_cuts = shortest_compatible_path_with_size_limit(
                    graph, 0, sequence_length, size_limit=path_size_limit,
                    min_step=path_size_min_step,
                    nodes_constraints=cuts_set_constraints
                )
            else:
                best_cuts = shortest_compatible_path(
                    graph, 0, sequence_length,
                    nodes_constraints=cuts_set_constraints
                )
        except (KeyError, nx.NetworkXNoPath):
            raise NoSolutionFoundError("Could not find a solution in "
                                       "optimize_cuts_with_graph")
    else:  # a_factor > 0 means use the A* algorithm

        graph.add_edges_from(segments)
        d = {}

        def compute_weight(start, end, props):
            """Compute the weight (cost) for segment (start, end).

            Parameter `props` is useless and is there for compatibility reasons
            """
            segment = tuple(sorted((start, end)))
            if segment in d:
                return d[segment]
            score = segment_score_function(segment)
            if score >= 0:
                result = score
            else:
                result = np.inf
            d[segment] = result
            return result

        try:
            best_cuts = astar_path(graph, 0, sequence_length,
                                   heuristic=lambda n1, n2: abs(
                                       a_star_factor * (n2 - n1)),
                                   weight=compute_weight)
        except (KeyError, nx.NetworkXNoPath):
            raise NoSolutionFoundError("Could not find a solution in "
                                       "optimize_cuts_with_graph")
    return graph, best_cuts


def refine_cuts_with_graph(sequence_length, cuts, radius,
                           segment_score_function, cut_location_constraints=(),
                           segment_constraints=(), nucleotide_resolution=1,
                           forced_cuts=(), progress_bars=True,
                           cuts_set_constraints=(),
                           a_star_factor=0):
    """Refines the cuts to optimize a cutting problem, using a local search.

    Given an initial list `cuts` of cutting locations, the method will try to
    replace each cut by a nearby location to see if it minimizes further the
    total segment score.

    Parameters
    ----------

    sequence_length
      Length of the Dna sequence to decompose

    cuts
      List of indices indicating the initial cutting pattern to improve

    radius
      Radius of the region to consider around each original cut for
      optimization

    segment_score_function
      Function `(start, end) -> score` yielding a score for a given sub-segment
      of the sequence (e.g. price of the segment)

    cut_location_constraints
      List or tuple of functions `int -> bool` which return for each location
      whether the location should be considered as a cutting site (True) or
      discarded (False) in the decomposition of the sequence.
      The locations considered are the locations which pass every filter
      in the `cut_location_constraints` list.

    segment_constraints
      List or tuple of functions `(int, int) -> bool` which return for each
      segment `(start, end)` whether the segment is a valid segment for the
      decomposition of the sequence or whether it should be forbidden.
      The segmentss considered are the segments which pass every filter
      in the `segments_filters` list.

    nucleotide_resolution
      Nucleotide resolution to use for the refinement. Not necessarily 1, but
      should be smaller than the nucleotide resolution used for the
      coarse-grain optimization.

    progress_bars
      If `progress_bars` is non-zero and `a_star_factor` is zero, a progress
      bar is shown as the edges of the problem are being evaluated

    a_star_factor
      If 0, the classical Dijkstra algorithm is used for path finding. Else,
      the A* algorithm is used with a heuristic `h(x)=a_start_factor*(L-x)`
      where x is the location of a cutting point and L the length of the
      sequence. See the original DNAWeaver article for clearer explanations.
      Using a high A* factor can improve computing times several folds but
      yields suboptimal decompositions.
      Note that using the A* algorithm doesnt play well with time constraints
      and number-of-fragments limits etc. Use with care

    Notes
    -----

    See original DNA Weaver article for clearer explanations of the method
    """
    nodes = [
        set([cut] + ([] if cut in forced_cuts else
                     list(range(max(0, cut - radius),
                                min(sequence_length + 1, cut + radius),
                                nucleotide_resolution))))
        for cut in cuts
    ]

    segments = [
        (start, end)
        for nodes_1, nodes_2 in zip(nodes, nodes[1:])
        for start, end in itt.product(nodes_1, nodes_2)
        if all(fl(start, end) for fl in segment_constraints)
    ]

    graph = nx.DiGraph()

    if a_star_factor == 0:
        for start, end in (segments if not progress_bars else
                           tqdm(segments, desc="Refining")):
            weight = segment_score_function((start, end))
            if weight >= 0:
                graph.add_edge(start, end, weight=weight)

        return shortest_compatible_path(
            graph, 0, sequence_length,
            nodes_constraints=cuts_set_constraints
        )
    else:
        graph.add_edges_from(segments)
        d = {}

        def compute_weight(start, end, props):
            segment = tuple(sorted((start, end)))
            if segment in d:
                return d[segment]
            score = segment_score_function(segment)
            if score >= 0:
                result = score
            else:
                result = np.inf
            d[segment] = result
            return result

        return astar_path(graph, 0, sequence_length,
                          heuristic=lambda n1, n2:
                              abs(a_star_factor * (n2 - n1)),
                          weight=compute_weight)


def optimize_cuts_with_graph_twostep(sequence_length,
                                     segment_score_function,
                                     cut_location_constraints=(),
                                     segment_constraints=(),
                                     initial_resolution=1,
                                     min_segment_length=500,
                                     max_segment_length=2000,
                                     refine_resolution=1,
                                     suggested_cuts=(),
                                     forced_cuts=(),
                                     cuts_set_constraints=(),
                                     a_star_factor=0,
                                     path_size_limit=None,
                                     path_size_min_step=0.001,
                                     progress_bars=True):
    """Find optimal sequence cuts in two steps: first coarse-grain search with
    a nucleotide resolution > 1, then a refinement step.

    This function calls `optimize_cuts_with_graph` followed by
    `refine_cuts_with_graph`.
    In practice, this means that at first only every N-th nucleotide of the
    sequence is considered as a possible cut point, then a refinement around
    the solution found is operated using local searches.


    Parameters
    ----------

    sequence_length
      Length of the sequence to decompose.

    segment_score_function
      Function `(start, end) -> score` yielding a score for a given sub-segment
      of the sequence (e.g. price of the segment)

    location_constraints
      List or tuple of functions `int -> bool` which return for each location
      whether the location should be considered as a cutting site (True) or
      discarded (False) in the decomposition of the sequence.
      The locations considered are the locations which pass every constraint
      in the `location_constraints` list.

    segment_constraints
      List or tuple of functions `(int, int) -> bool` which return for each
      segment `(start, end)` whether the segment is a valid segment for the
      decomposition of the sequence or whether it should be forbidden.
      The segmentss considered are the segments which pass every constraint
      in the `segments_constraints` list.

    initial_resolution
      Nucleotide resolution to use for the coarse-grain search, see function
      `optimize_cuts_with_graph`.

    min_segment_length
      Min length of the segments (this will be translated as an additional
      constraint in `segments_constraints`)

    max_segment_length
       Max length of the segments (this will be translated as an additional
       constraint in `segments_constraints`)

    refine_resolution
      nucleotide resolution to use during the local refinement step (see
      `refine_cuts_with_graph`).

    forced_cuts
      List of locations at which the decomposition must imperatively cut, even
      if these cuts do not comply with the `location_constraints`.

    a_star_factor
      If 0, the classical Dijkstra algorithm is used for path finding. Else,
      the A* algorithm is used with a heuristic `h(x)=a_start_factor*(L-x)`
      where x is the location of a cutting point and L the length of the
      sequence. See the original DNAWeaver article for clearer explanations.
      Using a high A* factor can improve computing times several folds but
      yields suboptimal decompositions.
      Note that using the A* algorithm doesnt play well with time constraints
      and number-of-fragments limits etc. Use with care

    path_size_limit
      Maximal number of edges for acceptable paths. None means no limit.
      Only works with `a_star_factor=0`

    path_size_min_step
      Minimal step for the bisection search when `path_size_limit` is not None
      (see `dijkstra_path_with_size_limit`)

    progress_bars
      If `progress_bars` is non-zero and `a_star_factor` is zero, a progress
      bar is shown as the edges of the problem are being evaluated

    Returns
    -------

    graph
      The graph of the problem (a Networkx DiGraph whose nodes are indices of
      locations in the sequence and the "weight" of edge [i][j] is given by
      segment_score_function(segment[i][j]).

    best_cuts
      The list of optimal cuts (includes 0 and len(sequence)), which minimizes
      the total score of all the segments corresponding to the cuts.

    Notes
    -----

    See original DnaWeaver article for more explanations

    """

    def is_resolution_location(location):
        """Return True iff the location is a N-th nucleotide, where N
        N is the nucleotide resolution of the problem."""
        return location % initial_resolution == 0

    new_location_constraints = ([is_resolution_location] +
                                list(cut_location_constraints))

    def size_is_valid(start, end):
        """Return True iff the segment's length is in interval
        [min_length, max_length]"""
        segment_length = end - start
        return min_segment_length <= segment_length <= max_segment_length

    new_segment_constraints = [size_is_valid] + list(segment_constraints)

    graph, best_cuts = optimize_cuts_with_graph(
        sequence_length,
        segment_score_function=segment_score_function,
        cut_location_constraints=new_location_constraints,
        segment_constraints=new_segment_constraints,
        forced_cuts=forced_cuts,
        cuts_set_constraints=cuts_set_constraints,
        suggested_cuts=suggested_cuts,
        max_segment_length=max_segment_length,
        a_star_factor=a_star_factor,
        progress_bars=progress_bars,
        path_size_limit=path_size_limit,
        path_size_min_step=path_size_min_step
    )
    if (initial_resolution > 1) and refine_resolution:
        best_cuts = refine_cuts_with_graph(
            sequence_length=sequence_length,
            cuts=best_cuts,
            radius=int(initial_resolution / 2),
            segment_score_function=segment_score_function,
            nucleotide_resolution=refine_resolution,
            segment_constraints=segment_constraints,
            cut_location_constraints=cut_location_constraints,
            forced_cuts=forced_cuts,
            cuts_set_constraints=cuts_set_constraints,
            a_star_factor=a_star_factor,
            progress_bars=progress_bars
        )
    return graph, best_cuts
