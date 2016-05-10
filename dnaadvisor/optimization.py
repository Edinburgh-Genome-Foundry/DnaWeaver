"""Optimization techniques"""

import numpy as np
from deap import creator, base, tools, algorithms
from tqdm import tqdm
import itertools as itt
import networkx as nx


def optimize_cuts_with_graph(sequence_length, segment_score_function,
                             cuts_number_penalty=0, location_filters=(),
                             segment_filters=()):

    nodes = sorted(list(set( [0, sequence_length] + [
        node
        for node in range(0, sequence_length)
        if all(fl(node) for fl in location_filters)
    ])))
    segments = [
        (start, end)
        for start, end in itt.product(nodes, nodes)
        if all(fl((start, end)) for fl in segment_filters)
    ]
    graph = nx.DiGraph()
    for start, end in tqdm(segments):
        weight = segment_score_function((start, end))
        if weight >= 0:
            graph.add_edge(start, end, weight=weight + cuts_number_penalty)
    best_cuts = nx.dijkstra_path(graph, 0, sequence_length)

    return best_cuts

def refine_cuts_with_graph(sequence_length, cuts, radius,
                           segment_score_function, location_filters=(),
                           segment_filters=(), nucleotide_resolution=1):
    nodes = [
        range(max(0, cut - radius),
              min(sequence_length + 1, cut + radius),
              nucleotide_resolution)
        for cut in cuts
    ]

    segments = [
        (start, end)
        for nodes_1, nodes_2 in zip(nodes, nodes[1:])
        for start, end in itt.product(nodes_1, nodes_2)
        if all(fl((start, end)) for fl in segment_filters)
    ]


    graph = nx.DiGraph()

    for start, end in tqdm(segments):
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
                                     refine_resolution=1):

    def is_resolution_location(location):
        return location % initial_resolution == 0

    new_location_filters = list(location_filters) + [is_resolution_location]

    def size_is_valid(segment):
        segment_length = segment[1] - segment[0]
        return min_segment_length < segment_length < max_segment_length

    new_segment_filters = list(segment_filters) + [size_is_valid]

    best_cuts = optimize_cuts_with_graph(
        sequence_length,
        segment_score_function=segment_score_function,
        cuts_number_penalty=cuts_number_penalty,
        location_filters=new_location_filters,
        segment_filters=new_segment_filters
    )
    if (initial_resolution > 1) and refine_resolution:
        best_cuts = refine_cuts_with_graph(
            sequence_length=sequence_length,
            cuts=best_cuts,
            radius=int(initial_resolution / 2),
            segment_score_function=segment_score_function,
            nucleotide_resolution=refine_resolution,
            segment_filters=segment_filters,
            location_filters=location_filters
        )
    return best_cuts
