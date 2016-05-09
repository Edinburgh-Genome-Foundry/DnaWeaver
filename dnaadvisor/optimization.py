"""Optimization techniques"""

import numpy as np
from deap import creator, base, tools, algorithms
from tqdm import tqdm
import itertools as itt
import networkx as nx


def optimize_cuts_with_graph(sequence, segment_score_function,
                             segment_length_range=(500, 2000),
                             nucleotide_resolution=10, refine=True):
    min_length, max_length = segment_length_range
    nodes = range(0, len(sequence), nucleotide_resolution)
    nodes[-1] = len(sequence)
    segments = [
        (start, end)
        for start, end in itt.product(nodes, nodes)
        if min_length < end - start < max_length
    ]
    graph = nx.DiGraph()
    for start, end in tqdm(segments):
        cost = segment_score_function((start, end), sequence)
        graph.add_edge(start, end, weight=cost)

    best_cuts = nx.dijkstra_path(graph, 0, len(sequence))

    if refine and (nucleotide_resolution > 1):
        radius = int(nucleotide_resolution / 2)
        best_cuts = refine_cuts_with_graph(sequence, best_cuts,
                                           segment_score_function, radius)
    return best_cuts

def refine_cuts_with_graph(sequence, cuts, segment_score_function, radius):
    nodes=[
        range(max(0, cut - radius), min(len(sequence)+1, cut + radius))
        for cut in cuts
    ]
    segments = [
        (node_1, node_2)
        for nodes_1, nodes_2 in zip(nodes, nodes[1:])
        for node_1, node_2 in itt.product(nodes_1, nodes_2)
        if node_2 - node_1 > 1
    ]
    graph=nx.DiGraph()
    for start, end in tqdm(segments):
        cost=segment_score_function((start, end), sequence)
        graph.add_edge(start, end, weight=cost)
    return nx.dijkstra_path(graph, 0, len(sequence))


def optimize_costs_with_graph(dna_ordering_problem,
                              segment_length_range=(500, 2000),
                              nucleotide_resolution=10, refine=True):
    def segment_score_function(segment, sequence):
        offers=dna_ordering_problem.find_offers(segment)
        if offers == []:
            return float(1e8)
        else:
            return min([offer.price for offer in offers])
    best_cuts = optimize_cuts_with_graph(
        dna_ordering_problem.sequence,
        segment_score_function,
        segment_length_range=segment_length_range,
        nucleotide_resolution=nucleotide_resolution,
        refine=refine
    )
    return dna_ordering_problem.find_best_offers(best_cuts)
