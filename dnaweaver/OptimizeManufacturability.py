from copy import deepcopy
import numpy as np
import dnachisel as dc

def compute_cost_profiles(graph):
    
    L = max(graph.nodes)
    cost_profiles = {
        'min': 1000 * np.ones(L),
        'max': -1000 * np.ones(L),
        'mean': np.zeros(L),
        'total': np.zeros(L),
        'count': np.zeros(L)
    }
    for start, end, data in graph.edges(data=True):
        cost_profiles['count'][start:end] += 1.0
        bp_price = data['weight'] / (end - start)
        cost_profiles['total'][start:end] += bp_price
        cost_profiles['min'][start:end] = np.minimum(
            bp_price, cost_profiles['min'][start:end])
        cost_profiles['max'][start:end] = np.maximum(
            bp_price, cost_profiles['max'][start:end])
    cost_profiles['mean'] = cost_profiles['total'] / (
        .0001 + cost_profiles['count'])
    diffs = np.diff(cost_profiles['min'])
    cost_profiles['diffs'] = diffs
    return cost_profiles

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
    
    def compute_cost_profiles(self, sequence):
        graph = self.compute_cost_graph(sequence)
        return compute_cost_profiles(graph)

    def detect_cost_driving_locations(self, sequence):
        diffs = self.compute_cost_profiles(sequence)['diffs']
        above_threshold = abs(diffs) > self.detection_threshold
        changes = [
            (index, diffs[index])
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