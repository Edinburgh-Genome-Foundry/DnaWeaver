from copy import deepcopy
import numpy as np
import dnachisel as dc

def compute_costs_profiles(graph):
    
    L = max(graph.nodes)
    costs_profiles = {
        'min': 1000 * np.ones(L),
        'max': -1000 * np.ones(L),
        'mean': np.zeros(L),
        'total': np.zeros(L),
        'count': np.zeros(L)
    }
    for start, end, data in graph.edges(data=True):
        costs_profiles['count'][start:end] += 1.0
        bp_price = data['weight'] / (end - start)
        costs_profiles['total'][start:end] += bp_price
        costs_profiles['min'][start:end] = np.minimum(
            bp_price, costs_profiles['min'][start:end])
        costs_profiles['max'][start:end] = np.maximum(
            bp_price, costs_profiles['max'][start:end])
    costs_profiles['mean'] = costs_profiles['total'] / (
        .0001 + costs_profiles['count'])
    diffs = np.diff(costs_profiles['min'])
    costs_profiles['diffs'] = diffs
    return costs_profiles

def plot_costs_profiles(bp_price_arrays):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    L = len(bp_price_arrays['count'])
    diffs = bp_price_arrays['diffs']
    r = (L - len(diffs)) / 2
    xx = range(L)
    ax1.plot(r + np.arange(len(diffs)), diffs, c='orange')
    ymax = diffs.max()
    if ymax == 0:
        ymax = 1
    ax1.set_ylim(bottom=-0.005, top = 1.5 * ymax)
    ax1b = ax1.twinx()
    ax1b.plot(xx, bp_price_arrays['count'], c='k', alpha=0.3)
    ax1b.set_ylim(bottom=0, top = 1.5 * bp_price_arrays['count'].max())
    ax1.set_xlim(0, L)
    
    ax2.plot(xx, bp_price_arrays['mean'], c='b', lw=2)
    ax2.fill_between(
        xx,
        bp_price_arrays['min'],
        bp_price_arrays['max'],
        facecolor='b',
        alpha=0.2
    )
    ax2.set_ylim(bottom=0, top = 1.5 * bp_price_arrays['max'].max())
    return ax1, ax2


class OptimizeManufacturability(dc.Specification):
    
    def __init__(self, station, detection_threshold=0.0005,
                 max_loc_size=50, resolution=10, boost=1):
        co_station = deepcopy(station)
        co_station.assembly_method.min_segment_length = resolution
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
        grained_cuts = set(range(0, decomposer.sequence_length + 1,
                           self.resolution))
        grained_cuts = grained_cuts.union(set([decomposer.sequence_length]))
        graph = decomposer.compute_graph(grained_cuts,
                                         reachable_indices_only=False)
        graph, cuts = decomposer.find_shortest_path(graph)
        return graph

    def detect_cost_driving_locations(self, sequence):
        graph = self.compute_cost_graph(sequence)
        costs_profiles = compute_costs_profiles(graph)
        threshold = self.detection_threshold
        above_threshold = abs(costs_profiles['diffs']) > threshold
        changes = [
            (index, costs_profiles['diffs'][index])
            for index in above_threshold.nonzero()[0]
        ]
        return [
            dc.Location(start[0], end[0])
            for start, end in zip(changes, changes[1:])
            if start[1] > 0
            and end[1] < 0
            and end[0] - start[0] < self.max_loc_size
        ]
        
    def evaluate(self, problem):
        if not self.is_localized:
            locations = self.detect_cost_driving_locations(problem.sequence)
        else:
            locations = []
        quote = self.get_quote(problem)
        score = -quote.price
        return dc.SpecEvaluation(self, problem=problem, score=score,
                                 locations=locations, data=quote)

    def get_quote(self, problem):
        return self.station.get_quote(problem.sequence,
                                      with_assembly_plan=False)

    def localized(self, location, problem=None):
        return self.copy_with_changes(is_localized=True)