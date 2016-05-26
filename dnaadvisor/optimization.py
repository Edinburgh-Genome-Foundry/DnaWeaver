"""Optimization techniques"""

import numpy as np
from deap import creator, base, tools, algorithms
from tqdm import tqdm
import itertools as itt
import networkx as nx


class NoSolutionFoundError(Exception):
    pass


def optimize_cuts_with_graph(sequence_length, segment_score_function,
                             cuts_number_penalty=0, location_filters=(),
                             segment_filters=(), forced_cuts=(),
                             max_segment_length=None, progress_bars=False):
    """Find the sequence cuts which optimize the sum of segments scores.

    This is a very generic method meant to be applied to any sequence cutting
    problem.

    Parameters
    ----------

    sequence_length
      Length of the sequence

    segment_score_function
      A function f( (start, end) ) -> score where (start, end) refers to the
      sequence segment start:end, and score is a float. The algorithm will
      produce cuts which minimize the total scores of the segments.

    cuts_number_penalty
      A penalty that can be applied


    location_filters

    segment_filters


    Returns
    -------

    graph
      The graph of the problem (a Networkx DiGraph whose nodes are indices of
      locations in the sequence and the "weight" of edge [i][j] is given by
      segment_score_function(segment[i][j] + cuts_number_penalty)).

    best_cuts
      The list of optimal cuts (includes 0 and len(sequence)), which minimizes
      the total score of all the segments corresponding to the cuts.

    """

    forced_cuts = list(forced_cuts)
    if max_segment_length is None:
        max_segment_length = sequence_length

    nodes = sorted(list(set([0, sequence_length] + [
        node
        for node in range(0, sequence_length)
        if all(fl(node) for fl in location_filters)
    ] + forced_cuts)))
    if forced_cuts != []:
        def forced_cuts_filter(segment):
            start, end = segment
            return not any((start < cut < end) for cut in forced_cuts)
        segment_filters = [forced_cuts_filter] + list(segment_filters)

    segments = []
    for i, start in enumerate(nodes if not progress_bars else
                              tqdm(nodes, desc="Filtering edges")):
        for end in nodes[i + 1:]:
            if end - start > max_segment_length:
                break
            elif all(fl((start, end)) for fl in segment_filters):
                segments.append((start, end))
    graph = nx.DiGraph()
    for start, end in (segments if not progress_bars else
                       tqdm(segments, desc='Computing edges')):
        weight = segment_score_function((start, end))
        if weight >= 0:
            graph.add_edge(start, end, weight=weight + cuts_number_penalty)

    try:
        best_cuts = nx.dijkstra_path(graph, 0, sequence_length)
    except (KeyError, nx.NetworkXNoPath) as err:
        raise NoSolutionFoundError("Could not find a solution in "
                                   "optimize_cuts_with_graph")

    return graph, best_cuts


def refine_cuts_with_graph(sequence_length, cuts, radius,
                           segment_score_function, location_filters=(),
                           segment_filters=(), nucleotide_resolution=1,
                           forced_cuts=(), progress_bars=True):

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
        if all(fl((start, end)) for fl in segment_filters)
    ]

    graph = nx.DiGraph()

    for start, end in (segments if not progress_bars else
                       tqdm(segments, desc="Refining")):
        weight = segment_score_function((start, end))
        if weight >= 0:
            graph.add_edge(start, end, weight=weight)

    return nx.dijkstra_path(graph, 0, sequence_length)


def optimize_cuts_with_graph_twostep(sequence_length,
                                     segment_score_function,
                                     cuts_number_penalty=0,
                                     location_filters=(),
                                     segment_filters=(),
                                     initial_resolution=1,
                                     min_segment_length=500,
                                     max_segment_length=2000,
                                     refine_resolution=1,
                                     forced_cuts=(),
                                     progress_bars=True):
    """Find optimal sequence cuts with coarse-grain search + refinement step.


    Parameters
    ----------

    initial_resolution

    min_segment_length

    max_segment_length

    refine_resolution

    other parameters
      Other parameters will be passed to functions optimize_cuts_with_graph and
      refine_cuts_with_graph.

    Returns
    -------

    graph
      The graph of the problem (a Networkx DiGraph whose nodes are indices of
      locations in the sequence and the "weight" of edge [i][j] is given by
      segment_score_function(segment[i][j] + cuts_number_penalty)).

    best_cuts
      The list of optimal cuts (includes 0 and len(sequence)), which minimizes
      the total score of all the segments corresponding to the cuts.

    """
    def is_resolution_location(location):
        return location % initial_resolution == 0

    new_location_filters = [is_resolution_location] + list(location_filters)

    def size_is_valid(segment):
        segment_length = segment[1] - segment[0]
        return min_segment_length < segment_length < max_segment_length

    new_segment_filters = [size_is_valid] + list(segment_filters)

    graph, best_cuts = optimize_cuts_with_graph(
        sequence_length,
        segment_score_function=segment_score_function,
        cuts_number_penalty=cuts_number_penalty,
        location_filters=new_location_filters,
        segment_filters=new_segment_filters,
        forced_cuts=forced_cuts,
        max_segment_length=max_segment_length,
        progress_bars=progress_bars
    )
    if (initial_resolution > 1) and refine_resolution:
        best_cuts = refine_cuts_with_graph(
            sequence_length=sequence_length,
            cuts=best_cuts,
            radius=int(initial_resolution / 2),
            segment_score_function=segment_score_function,
            nucleotide_resolution=refine_resolution,
            segment_filters=segment_filters,
            location_filters=location_filters,
            forced_cuts=forced_cuts,
            progress_bars=progress_bars
        )
    return graph, best_cuts
