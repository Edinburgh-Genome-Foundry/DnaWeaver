import numpy as np
from functools import lru_cache
from .SegmentSelector import SegmentSelector

try:
    import primer3

    PRIMER3_AVAILABLE = True
except ImportError:
    PRIMER3_AVAILABLE = False


class TmSegmentSelector(SegmentSelector):
    """Selects segment with melting temperature constraints.

    Given a location in a sequence, this selector will return a segment
    (=subsegment of the sequence) around this position, with size and
    temperature constraints.

    Can be used for Gibson Assembly, Golden Gate Assembly, etc.

    This class also implements methods to quickly compute all segments in a
    sequence. The result is then cached for speed.

    Parameters
    ----------

    min_size
      Minimal length of the segment, in nucleotides.

    max_size
      Maximal length of the segment, in nucleotides.

    min_tm
      Minimal melting temp allowed for the segments.

    max_tm
      Maximal melting temp allowed for the segments.

    precompute_segments
      If True, the whole sequence will be analyzed to find the best segment
      for each cutting point (and find cutting points which do not allow
      segments), using the AT/GC=2/4C heuristic for melting temperatures.

    primer3_params
      If ``precompute_segments`` is False and ``Primer3`` is installed, the
      melting temperature is computed.
    """

    def __init__(
        self,
        min_size=18,
        max_size=22,
        min_tm=55,
        max_tm=65,
        precompute_segments=True,
        left_addition="",
        right_addition="",
        **primer3_params
    ):
        """Initialize."""
        self.min_size = min_size
        self.max_size = max_size
        self.precompute_segments = precompute_segments
        self.primer3_params = primer3_params
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.middle_tm = 0.5 * (self.min_tm + self.max_tm)
        self.middle_size = int(0.5 * (self.min_size + self.max_size))
        self.left_addition = left_addition
        self.right_addition = right_addition

    @property
    def max_homology_size(self):
        return self.max_size

    def filter_location(self, sequence, index):
        """Return whether the sequence has a valid segment at this index."""
        if index == len(sequence):
            index = index - 1
        if self.precompute_segments:
            return self.compute_all_segments(sequence)[index] is not None
        else:
            return self.compute_optimal_segment(sequence, index)[1]

    def compute_segment_location(self, sequence, index):
        """Return the location (start, stop) of the selected segment."""
        if (index < 0) or (sequence == ""):
            return None
        if index == len(sequence):
            index = index - 1
        if self.precompute_segments:
            # calling a result that should be computed once then cached
            return self.compute_all_segments(sequence)[index]
        else:
            return self.compute_optimal_segment(sequence, index)[0]

    def compute_optimal_segment(self, sequence, index, init_size=None):
        """Return the best segment around some location.

        This function will be shortcut if ``precompute_segments`` is True,
        as it is slow to compute the best overhang for each location
        separately.
        """

        if init_size is None:
            init_size = self.middle_size

        def f(ovh_size):
            start, end = self.get_segment_coordinates(index, ovh_size, len(sequence))
            tm = self.compute_tm(sequence[start:end])
            return (start, end), tm

        best_segment, best_tm = f(init_size)
        if best_tm < self.middle_tm:
            ovh_sizes = range(init_size, self.max_size)
        else:
            ovh_sizes = range(init_size, self.min_size, -1)
        for ovh_size in ovh_sizes:
            new_segment, new_tm = f(ovh_size)
            if abs(new_tm - self.middle_tm) < abs(best_tm - self.middle_tm):
                best_segment, best_tm = new_segment, new_tm
            else:
                break
        return best_segment, (self.min_tm < best_tm < self.max_tm)

    @lru_cache(maxsize=3)
    def compute_all_segments(self, sequence):
        """Quickly compute all segments in the sequence.

        This function uses the heuristic {A, T}=2degC, {G, C}=4degC to compute
        melting temperatures.

        This function uses vectorial operations for speed. The results are also
        cached.
        """

        if self.primer3_params == {}:
            lmin, lmax = self.min_size, self.max_size
            table = np.zeros((lmax + 1 - lmin, len(sequence)))

            # super fast GC computing using the trick 67=C, 71=G in unicode.
            seq_array = np.frombuffer((sequence + "").encode(), dtype="uint8")
            cumsum = np.cumsum(2 + 2 * ((seq_array == 71) | (seq_array == 67)))
            for i, oh_size in enumerate(range(lmin, min(len(sequence) - 1, lmax + 1))):
                arr = cumsum[oh_size:] - cumsum[:-oh_size]
                start = int(oh_size / 2)
                end = start + len(arr)
                table[i, start:end] = arr
                # print (i, start, table)
                table[i, :start] = table[i, start]
                table[i, end:] = table[i, end - 1]
            scores = -(table - self.min_tm) * (table - self.max_tm)
            best_sizes_indices = scores.argmax(axis=0)
            best_sizes = lmin + best_sizes_indices
            validities = np.choose(best_sizes_indices, scores) >= 0
            osizes_and_validities = zip(best_sizes, validities)
        else:
            osizes_and_validity = [self.compute_optimal_segment(sequence, 0)]
            for i in range(1, len(sequence)):
                init_size = osizes_and_validity[-1][0]
                ovh, val = self.find_optimal_segment(sequence, i, init_size)
                osizes_and_validity.append((len(ovh), val))
        return [
            None
            if not valid
            else self.get_segment_coordinates(i, ovh_size, len(sequence))
            for i, (ovh_size, valid) in enumerate(osizes_and_validities)
        ]

    def compute_tm(self, sequence):
        """Return the melting temp of the sequence.

        If Primer3 is available, it's internal melting temperature calculator
        is used with ``self.primer3_params`` used as parameters.
        Else the heuristic AT/GC=2/4C is used.
        """
        if self.primer3_params == {}:
            return sum([4 if c in "GC" else 2 for c in sequence])
        if not PRIMER3_AVAILABLE:
            raise ImportError(
                "Melting temperature computation with '%s' "
                "Requires Primer3 installed."
                % self.primer3_params.get("method", "[unknown method]")
            )
        return primer3.calcTm(sequence, **self.primer3_params)

    def __str__(self):
        result = "Tm(%d-%dC, %d-%dbp)" % (
            self.min_tm,
            self.max_tm,
            self.min_size,
            self.max_size,
        )
        if self.left_addition:
            result = ("...%s-" % self.left_addition[-12:]) + result
        if self.right_addition:
            result = result + ("-%s..." % self.right_addition[:12])
        return result
