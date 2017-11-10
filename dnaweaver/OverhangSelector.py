""" This module contains overhangs selectors for methods such as
Gibson Assembly, Golden Gate assembly, recombination in yeast."""

import numpy as np
from functools import lru_cache
try:
    import primer3
    PRIMER3_AVAILABLE = True
except:
    PRIMER3_AVAILABLE = False



class OverhangSelector:
    """Base class for overhang selectors such as TmOverhangSelector.

    These selectors return an ideal overhang subsegment around a given sequence
    location. They can also be used to filter out locations in the sequence as
    non-potential cutting sites for sequence decomposition.
    """

    def location_filter_method(self, sequence):
        """Return a filter function f(location) => True/False.

        The result is True iff the location is a valid cutting site"""
        def f(index):
            return self.filter_location(sequence, index)
        return f

    def filter_location(self, sequence, index):
        """By default, an OverhangSelector does not filter out any location."""
        return True

    def compute_sequence_fragment(self, sequence, segment):
        """Compute the fragment equal to the sequence segment + overhangs."""
        start, end = segment
        if start == 0:
            fragment_start = 0
        else:
            fragment_start = self.compute_overhang_location(sequence, start)[0]
        if end == len(sequence):
            fragment_end = len(sequence)
        else:
            fragment_end = self.compute_overhang_location(sequence, end)[1]
        fragment = sequence[fragment_start: fragment_end]
        return self.left_addition + fragment + self.right_addition


    def compute_overhang_sequence(self, sequence, index):
        """Return the sequence of the selected overhang at the given index."""
        start, end = self.compute_overhang_location(sequence, index)
        return sequence[start: end]

    @staticmethod
    def get_segment_coordinates(center, segment_length, sequence_length):
        """Return max(0, c - s/2) - min(L, c + L/2).

        Where c=center, s=segment_length, L=sequence_length.
        """
        half = int(segment_length / 2)
        start = max(0, min(center - half, sequence_length - segment_length))
        end = start + segment_length
        return start, end


class FixedSizeOverhangSelector(OverhangSelector):
    """Selects overhangs of a constant size.

    Great for methods involving large homology regions where melting
    temperature matters less.
    """

    def __init__(self, overhang_size=100, left_addition='', right_addition=''):
        self.overhang_size = overhang_size
        self.left_addition = left_addition
        self.right_addition = right_addition


    def compute_overhang_location(self, sequence, index):
        return self.get_segment_coordinates(index, self.overhang_size,
                                            len(sequence))

class TmOverhangSelector(OverhangSelector):
    """Selects overhangs with melting temperature constraints.

    Given a location in a sequence, this selector will return an overhang
    (subsegment of the sequence) around this position, with size and
    temperature constraints.

    Can be used for Gibson Assembly, Golden Gate Assembly, etc.

    This class also implements methods to quickly compute all overhangs in a
    sequence. The result is then cached for speed.

    Parameters
    ----------
    min_size
      Minimal length of the overhang, in nucleotides

    max_size
      Maximal length of the overhang, in nucleotides

    min_tm
      Minimal melting temp allowed for overhangs

    max_tm
      Maximal melting temp allowed for overhangs

    precompute_overhangs
      If True, the whole sequence will be analyzed to find the best overhang
      for each cutting point (and find cutting points which do not allow
      overhangs), using the AT/GC=2/4C heuristic for melting temperatures.

    primer3_params
      If ``precompute_overhangs`` is False and ``Primer3`` is installed, the
      melting temperature is computed
    """

    def __init__(self, min_size=18, max_size=22, min_tm=55, max_tm=65,
                 precompute_overhangs=True, left_addition='',
                 right_addition='', **primer3_params):
        """Initialize."""
        self.min_size = min_size
        self.max_size = max_size
        self.precompute_overhangs = precompute_overhangs
        self.primer3_params = primer3_params
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.middle_tm = 0.5 * (self.min_tm + self.max_tm)
        self.middle_size = int(0.5 * (self.min_size + self.max_size))
        self.left_addition = left_addition
        self.right_addition = right_addition

    def filter_location(self, sequence, index):
        """Return whether the sequence has a valid overhang at this index."""
        if index == len(sequence):
            index = index - 1
        if self.precompute_overhangs:
            return (self.compute_all_overhangs(sequence)[index] is not None)
        else:
            return self.compute_optimal_overhang(sequence, index)[1]

    def compute_overhang_location(self, sequence, index):
        """Return the location (start, stop) of the selected overhang."""
        if index == len(sequence):
            index = index - 1
        if self.precompute_overhangs:
            # calling a result that should be computed once then cached
            return self.compute_all_overhangs(sequence)[index]
        else:
            return self.compute_optimal_overhang(sequence, index)[0]

    def compute_optimal_overhang(self, sequence, index, init_size=None):
        """Return the best overhang around some location.

        This function will be shortcut if ``precompute_overhangs`` is True,
        as it is slow to compute the best overhang for each location
        separately.
        """

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

    @lru_cache(maxsize=3)
    def compute_all_overhangs(self, sequence):
        """Quickly compute all overhangs in the sequence.

        This function uses the heuristic {A, T}=2degC, {G, C}=4degC to compute
        melting temperatures.

        This function uses vectorial operations for speed. The results are also
        cached.
        """

        if self.primer3_params == {}:
            lmin, lmax = self.min_size, self.max_size
            table = np.zeros((lmax + 1 - lmin, len(sequence)))

            # super fast GC computing using the trick 67=C, 71=G in unicode.
            seq_array = np.fromstring(sequence + "", dtype="uint8")
            cumsum = np.cumsum(2 + 2 * ((seq_array == 71) | (seq_array == 67)))
            # cumsum = np.cumsum([
            #     4 if nuc in "GC" else 2
            #     for nuc in sequence
            # ])
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
        """Return the melting temp of the sequence.

        If Primer3 is available, it's internal melting temperature calculator
        is used with ``self.primer3_params`` used as parameters.
        Else the heuristic AT/GC=2/4C is used.
        """
        if self.params == {}:
            return sum([
                4 if c in "GC" else 2
                for c in sequence
            ])
        if not PRIMER3_AVAILABLE:
            raise ImportError("Melting temperature computation with '%s' "
                              "Requires Primer3 installed." %
                              self.primer3_params.get('method',
                                                      "[unknown method]"))
        return primer3.calcTm(sequence, **self.primer3_params)
