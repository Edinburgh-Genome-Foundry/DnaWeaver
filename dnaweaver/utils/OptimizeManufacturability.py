"""Module to optimize sequences for manufacturability using DnaWeaver.

This is still very experimental!
"""

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import dnachisel as dc
from dnachisel.biotools import group_nearby_indices
from ..DnaAssemblyMethod import BluntEndAssemblyMethod


class OptimizeManufacturability(dc.Specification):
    def __init__(
        self,
        station,
        detection_threshold=0.0001,
        max_loc_size=50,
        resolution=10,
        max_segment_length=5000,
        boost=1,
    ):
        co_station = deepcopy(station)
        co_station.assembly_method = BluntEndAssemblyMethod(
            min_segment_length=5 * resolution, max_segment_length=max_segment_length,
        )
        co_station.a_star_factor = 0
        self.station = station
        self.co_station = co_station
        self.resolution = resolution
        self.detection_threshold = detection_threshold
        self.boost = boost
        self.max_loc_size = max_loc_size
        self.is_localized = False

    def compute_cost_graph(self, sequence):
        decomposer = self.co_station.new_sequence_decomposer(sequence)
        grained_cuts = set(range(0, decomposer.sequence_length + 1, self.resolution))
        grained_cuts = grained_cuts.union(set([decomposer.sequence_length]))
        graph = decomposer.compute_graph(grained_cuts, reachable_indices_only=False)
        graph, _ = decomposer.find_shortest_path(graph)
        return graph

    def compute_cost_profiles(self, sequence):
        graph = self.compute_cost_graph(sequence)
        return compute_cost_profiles(graph)

    def detect_cost_driving_locations(self, sequence):
        costs = self.compute_cost_profiles(sequence)
        diffs = np.diff(costs["mean"][:: self.resolution], 2)
        locations = (diffs < -self.detection_threshold).nonzero()[0]
        locations = (locations + 1) * self.resolution
        locations = group_nearby_indices(
            locations, max_gap=self.max_loc_size, max_group_spread=self.max_loc_size,
        )
        return [
            dc.Location(l[0] - self.resolution, l[-1] + self.resolution)
            for l in locations
        ]

    def evaluate(self, problem):
        if not self.is_localized:
            locations = self.detect_cost_driving_locations(problem.sequence)
        else:
            locations = []
        quote = self.get_quote(problem)
        if not quote.accepted:
            score = -1e9
        else:
            score = -quote.price
        return dc.SpecEvaluation(
            self, problem=problem, score=score, locations=locations, data=quote
        )

    def get_quote(self, problem):
        return self.station.get_quote(problem.sequence, with_assembly_plan=False)

    def localized(self, location, problem=None):
        return self.copy_with_changes(is_localized=True)

    @staticmethod
    def plot_cost_profiles(cost_profiles):
        """Plot cost profiles generated in OptimizeManufacturability"""
        _fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        L = len(cost_profiles["count"])
        diffs = cost_profiles["diffs"]
        r = (L - len(diffs)) / 2
        xx = range(L)
        ax1.plot(r + np.arange(len(diffs)), diffs, c="orange")
        ymax = diffs.max()
        if ymax == 0:
            ymax = 1
        ax1.set_ylim(bottom=-0.005, top=1.5 * ymax)
        ax1b = ax1.twinx()
        ax1b.plot(xx, cost_profiles["count"], c="k", alpha=0.3)
        ax1b.set_ylim(bottom=0, top=1.5 * cost_profiles["count"].max())
        ax1.set_xlim(0, L)
        ax1.set_title("Basepair cost jumps")

        ax2.plot(xx, cost_profiles["mean"], c="b", lw=2)
        ax2.fill_between(
            xx, cost_profiles["min"], cost_profiles["max"], facecolor="b", alpha=0.2,
        )
        ax2.set_ylim(bottom=0, top=1.5 * cost_profiles["max"].max())
        ax2.set_title("<min cost, AVERAGE, max cost>")
        return ax1, ax2


def compute_cost_profiles(graph):

    L = max(graph.nodes)
    cost_profiles = {
        "min": 1000 * np.ones(L),
        "max": -1000 * np.ones(L),
        "mean": np.zeros(L),
        "total": np.zeros(L),
        "count": np.zeros(L),
    }
    for start, end, data in graph.edges(data=True):
        cost_profiles["count"][start:end] += 1.0
        bp_price = data["weight"] / (end - start)
        cost_profiles["total"][start:end] += bp_price
        cost_profiles["min"][start:end] = np.minimum(
            bp_price, cost_profiles["min"][start:end]
        )
        cost_profiles["max"][start:end] = np.maximum(
            bp_price, cost_profiles["max"][start:end]
        )
    cost_profiles["mean"] = cost_profiles["total"] / (0.0001 + cost_profiles["count"])
    diffs = np.diff(cost_profiles["mean"])
    cost_profiles["diffs"] = diffs
    return cost_profiles
