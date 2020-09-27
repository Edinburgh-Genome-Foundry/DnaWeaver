"""Optimization techniques."""

import numpy as np
import networkx as nx
from .shortest_path_algorithms import (
    NoSolutionFoundError,
    shortest_valid_path,
    astar_path,
)
from proglog import default_bar_logger

# _default_bars = ('segment', 'edge')
# class SequenceDecomposerLogger(TqdmProgressBarLogger):

#     def __init__(self, bars=_default_bars, notebook='default',
#                  min_time_interval=0.2):
#         ignored_bars = set(_default_bars).difference(bars)
#         TqdmProgressBarLogger.__init__(self, bars=bars, notebook=notebook,
#                                        ignored_bars=ignored_bars,
#                                        min_time_interval=min_time_interval)


class SequenceDecomposer:
    """Find the sequence cuts which optimize the sum of segment scores.

    This is a very generic method meant to be applied to any sequence
    cutting optimization problem.

    Parameters
    ----------

    sequence_length
      Length of the sequence.

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
      Only works with `a_star_factor=0`.

    path_size_min_step
      Minimal step for the bisection search when `path_size_limit` is not None
      (see `dijkstra_path_with_size_limit`).


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

    def __init__(
        self,
        sequence_length,
        segment_score_function,
        cut_location_constraints=(),
        segment_constraints=(),
        forced_cuts=(),
        suggested_cuts=(),
        cuts_set_constraints=(),
        cut_spread_radius=0,
        min_segment_length=0,
        max_segment_length=None,
        coarse_grain=1,
        fine_grain=1,
        a_star_factor=0,
        path_size_limit=None,
        path_size_min_step=0.001,
        logger=None,
        bar_prefix="",
    ):

        self.segment_score_function = segment_score_function
        self.forced_cuts = set(list(forced_cuts) + [0, sequence_length])
        self.suggested_cuts = set(suggested_cuts)
        self.cut_spread_radius = cut_spread_radius
        self.segments_constraints = list(segment_constraints)
        self.cuts_set_constraints = list(cuts_set_constraints)
        self.cut_location_constraints = list(cut_location_constraints)
        self.sequence_length = sequence_length
        if max_segment_length is None:
            self.max_segment_length = sequence_length
        else:
            self.max_segment_length = max_segment_length
        self.min_segment_length = min_segment_length

        self.coarse_grain = coarse_grain
        self.fine_grain = fine_grain
        self.a_star_factor = a_star_factor
        self.path_size_limit = path_size_limit
        self.path_size_min_step = path_size_min_step
        self.logger = default_bar_logger(logger, min_time_interval=0.2)
        self.bar_prefix = bar_prefix

        if len(forced_cuts) > 0:

            def forced_cuts_filter(segment):
                start, end = segment
                return not any((start < cut < end) for cut in forced_cuts)

            self.segments_constraints.append(forced_cuts_filter)

        if len(self.cut_location_constraints) == 0:
            self.valid_cuts = set(range(sequence_length + 1))
        else:
            self.valid_cuts = set(
                [
                    index
                    for index in range(sequence_length + 1)
                    if all(fl(index) for fl in self.cut_location_constraints)
                ]
            )

    def compute_graph(
        self, valid_cuts, prune_deadends=True, reachable_indices_only=True
    ):
        L = self.sequence_length
        segments = []
        reachable_indices = set({0})

        # LIST THE SEGMENTS
        valid_cuts = sorted(valid_cuts)
        valid_cuts_array = np.zeros(L + 1)
        valid_cuts_array[valid_cuts] = 1
        for start in self.logger.iter_bar(
            segment=valid_cuts, bar_prefix=self.bar_prefix
        ):
            if reachable_indices_only and (start not in reachable_indices):
                continue
            ends_min = start + self.min_segment_length
            ends_max = min(L + 1, start + self.max_segment_length)
            hits = valid_cuts_array[ends_min:ends_max].nonzero()[0]
            valid_ends = [int(e) for e in hits + ends_min]
            if len(self.segments_constraints) > 0:
                valid_ends = [
                    end
                    for end in valid_ends
                    if all(fl((start, end)) for fl in self.segments_constraints)
                ]
            segments += [(start, end) for end in valid_ends]
            reachable_indices = reachable_indices.union(valid_ends)
        graph = nx.DiGraph(segments)

        if prune_deadends:
            ancestors = nx.ancestors(graph, L)
            graph.remove_nodes_from(
                [n for n in graph if (n != L) and (n not in ancestors)]
            )
        return graph

    def find_shortest_path(self, graph):
        constraints = self.cuts_set_constraints
        if (len(constraints) == 0) and (self.a_star_factor > 0):
            memodict = {}

            def compute_weight(start, end, props):
                """Compute the weight (cost) for segment (start, end).

                Parameter `props` is useless and is there for compatibility
                reasons.
                """
                segment = tuple(sorted((start, end)))
                if segment in memodict:
                    return memodict[segment]
                score = self.segment_score_function(segment)
                if score >= 0:
                    result = score
                else:
                    result = np.inf
                memodict[segment] = result
                return result

            def heuristic(n1, n2):
                return self.a_star_factor * (n2 - n1)

            try:
                path = astar_path(
                    graph,
                    0,
                    self.sequence_length,
                    heuristic=heuristic,
                    weight=compute_weight,
                )
                return graph, path
            except (KeyError, nx.NetworkXNoPath):
                return graph, None

        else:
            for start, end in self.logger.iter_bar(
                edge=list(graph.edges()), bar_prefix=self.bar_prefix
            ):
                weight = self.segment_score_function((start, end))
                if weight < 0:
                    graph.remove_edge(start, end)
                else:
                    graph[start][end]["weight"] = weight

            try:
                path = shortest_valid_path(
                    graph,
                    0,
                    self.sequence_length,
                    nodes_constraints=constraints,
                    size_limit=self.path_size_limit,
                    min_step=self.path_size_min_step,
                )
                return graph, path
            except (KeyError, nx.NetworkXNoPath):
                return graph, None

    def compute_coarse_cuts(self, grain="default"):
        L = self.sequence_length
        if grain == "default":
            grain = self.coarse_grain
        grained_cuts = set(range(0, L + 1, grain))
        forced_but_invalid = self.forced_cuts.difference(self.valid_cuts)
        if len(forced_but_invalid):
            raise NoSolutionFoundError(
                "One forced cut was also invalid at indices %s" % forced_but_invalid
            )
        grained_and_suggested = grained_cuts.union(self.suggested_cuts)
        if self.cut_spread_radius:
            grained_and_suggested = spread_indices_in_set(
                indices=grained_and_suggested,
                radius=self.cut_spread_radius,
                max_index=L + 1,
            )
        valid_cuts = grained_and_suggested.intersection(self.valid_cuts).union(
            self.forced_cuts
        )

        self.coarse_graph = self.compute_graph(valid_cuts)
        error = NoSolutionFoundError(
            "No coarse solution found, possibly too "
            "strong cuts/segments constraints."
        )
        if 0 not in self.coarse_graph.nodes():
            raise error
        graph, cuts = self.find_shortest_path(self.coarse_graph)

        self.weighted_coarse_graph = graph
        if cuts is None:
            raise error
        self.coarse_cuts = cuts

    def compute_fine_cuts(self):
        L = self.sequence_length
        radius = int(self.coarse_grain / 2)
        fine_cuts = set().union(
            *[
                set(
                    [cut]
                    + (
                        []
                        if cut in self.forced_cuts
                        else list(
                            range(
                                max(0, cut - radius),
                                min(L + 1, cut + radius),
                                self.fine_grain,
                            )
                        )
                    )
                )
                for cut in self.coarse_cuts
            ]
        )
        if self.cut_spread_radius:
            fine_cuts = spread_indices_in_set(
                indices=fine_cuts, radius=self.cut_spread_radius, max_index=L + 1,
            )
        fine_cuts = (fine_cuts.intersection(self.valid_cuts)).union(self.forced_cuts)
        self.fine_graph = self.compute_graph(fine_cuts)
        _, self.fine_cuts = self.find_shortest_path(self.fine_graph)

    def compute_optimal_cuts(self):
        self.compute_coarse_cuts()
        if (self.coarse_grain > 1) and self.fine_grain:
            self.compute_fine_cuts()
            return self.fine_cuts
        else:
            return self.coarse_cuts


def spread_indices_in_set(indices, radius, max_index):
    new_list = [
        new_index
        for i in indices
        for new_index in range(max(0, i - radius), min(max_index, i + radius + 1))
    ]
    return set(new_list)
