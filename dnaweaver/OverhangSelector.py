import numpy as np
from cachetools import LRUCache, cached
try:
    import primer3
    PRIMER3_AVAILABLE = True
except:
    PRIMER3_AVAILABLE = False



class OverhangSelector:

    def location_filter_method(self, sequence):
        def f(index):
            return self.filter_location(sequence, index)
        return f

    def filter_location(self, sequence, index):
        return True

    def compute_sequence_fragment(self, sequence, segment):
        start, end = segment
        if start == 0:
            fragment_start = 0
        else:
            fragment_start = self.compute_overhang_location(sequence, start)[0]
        if end == len(sequence):
            fragment_end = len(sequence)
        else:
            fragment_end = self.compute_overhang_location(sequence, end)[1]
        return sequence[fragment_start: fragment_end]


    def compute_overhang_sequence(self, sequence, index):
        start, end = self.compute_overhang_location(sequence, index)
        return sequence[start: end]

    @staticmethod
    def get_segment_coordinates(center, segment_length, sequence_length):
        half = int(segment_length/2)
        start = max(0, min(center - half, sequence_length - segment_length))
        end = start + segment_length
        return start, end


class ConstantSizeOverhangSelector(OverhangSelector):

    def __init__(self, overhang_size=100):
        self.overhang_size = overhang_size

    def compute_overhang_location(self, sequence, index):
        return self.get_segment_coordinates(index, self.overhang_size,
                                            len(sequence))

class TmOverhangSelector(OverhangSelector):

    def __init__(self, min_size=18, max_size=22, min_tm=55, max_tm=65,
                 precompute_overhangs=True, **params):
        self.min_size = min_size
        self.max_size = max_size
        self.precompute_overhangs = precompute_overhangs
        self.params = params
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.middle_tm = 0.5 * (self.min_tm + self.max_tm)
        self.middle_size = int(0.5 * (self.min_size + self.max_size))

    def filter_location(self, sequence, index):
        if index == len(sequence):
            index = index - 1
        if self.precompute_overhangs:
            # calling a result that should be computed once then cached
            # bools = np.array([0.0 if e is None else 1.0
            #                   for e in self.compute_all_overhangs(sequence)])
            # print (bools.mean())

            return self.compute_all_overhangs(sequence)[index] is not None
        else:
            return self.compute_optimal_overhang(sequence, index)[1]

    def compute_overhang_location(self, sequence, index):
        if index == len(sequence):
            index = index - 1
        if self.precompute_overhangs:
            # calling a result that should be computed once then cached
            return self.compute_all_overhangs(sequence)[index]
        else:
            return self.compute_optimal_overhang(sequence, index)[0]

    def compute_optimal_overhang(self, sequence, index, init_size=None):
        if init_size is None:
            init_size = self.middle_size
        def f(ovh_size):
            start, end = self.get_segment_coordinates(index, ovh_size,
                                                      len(sequence))
            tm = self.compute_tm(sequence[start: end])
            return (start, end), tm
        best_overhang, best_tm = f(init_size)
        if best_tm < self.middle_tm:
            ovh_sizes = range(init_size, self.max_size)
        else:
            ovh_sizes = range(init_size, self.min_size, -1)
        for ovh_size in ovh_sizes:
            new_overhang, new_tm = f(ovh_size)
            if abs(new_tm - self.middle_tm) < abs(best_tm - self.middle_tm):
                best_overhang, best_tm = new_overhang, new_tm
            else:
                break
        return best_overhang, (self.min_tm < best_tm < self.max_tm)

    @cached(LRUCache(3))
    def compute_all_overhangs(self, sequence):
        if self.params == {}:
            lmin, lmax = self.min_size, self.max_size
            table = np.zeros((lmax + 1 - lmin, len(sequence)))
            cumsum = np.cumsum([
                2 if nuc in "GC" else 4
                for nuc in sequence
            ])
            for i, oh_size in enumerate(range(lmin, lmax + 1)):
                arr = cumsum[oh_size:] - cumsum[:-oh_size]
                start = int(oh_size / 2)
                end = start + len(arr)
                table[i, start:end] = arr
                table[i, :start] = table[i, start]
                table[i, end:] = table[i, end-1]
            scores = - (table - self.min_tm) * (table - self.max_tm)
            best_sizes_indices = scores.argmax(axis=0)
            best_sizes = lmin + best_sizes_indices
            validities = np.choose(best_sizes_indices, scores) >= 0
            osizes_and_validities = zip(best_sizes, validities)
        else:
            osizes_and_validity = [self.compute_optimal_overhang(sequence, 0)]
            for i in range(1, len(sequence)):
                init_size = osizes_and_validity[-1][0]
                ovh, val = self.find_optimal_overhang(sequence, i, init_size)
                osizes_and_validity.append((len(ovh), val))

        return [
          None if not valid
          else self.get_segment_coordinates(i, ovh_size, len(sequence))
          for i, (ovh_size, valid) in enumerate(osizes_and_validities)
        ]


    def compute_tm(self, sequence):
        if self.params == {}:
            return sum([
                4 if c in "GC" else 2
                for c in sequence
            ])
        if not PRIMER3_AVAILABLE:
            raise ImportError("Melting temperature computation with '%s' "
                              "Requires Primer3 installed." %
                              self.params.get('method', "[unknown method]"))
        return primer3.calcTm(sequence, **self.params)
